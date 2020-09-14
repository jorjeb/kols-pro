#!/usr/local/bin/python
from datetime import date
from igraph import Graph, operator, VertexSeq, WEAK
from flask import Flask, render_template, request, redirect
from neo4j.v1 import GraphDatabase, basic_auth
from wsgiref.handlers import CGIHandler
import csv
import gearman
import json
import os
import re
import sys
import uuid

reload(sys)
sys.setdefaultencoding('utf8')

app = Flask(__name__)
data_dir = '../data'  # no trailing slash
gearman_url = '127.0.0.1:4730'
neo4j_url = 'bolt://localhost'
neo4j_user = 'neo4j'
neo4j_password = '1234'


def clean_country(c):
    try:
        pos = c.rindex('(')
        return c[0:pos - 1]
    except ValueError:
        return c


def unescape(s):
    return re.sub(u'(^`|`$)', '', s)


def make_safe(f):
    f = f.replace(' ', '_')
    return "".join([c for c in f if re.match(r'[\w\-]', c)])


class DataEncoder(gearman.DataEncoder):
    @classmethod
    def encode(cls, object):
        return json.dumps(object)

    @classmethod
    def decode(cls, json_string):
        return json.loads(json_string)


class Client(gearman.GearmanClient):
    data_encoder = DataEncoder


@app.route('/')
def index():
    f = make_safe(request.args.get('f')) if request.args.get('f') else ''
    results = []

    if f:
        filename = "{}/{}.csv".format(data_dir, f)

        if os.path.isfile(filename):
            with open(filename, 'r') as f:
                reader = csv.reader(f, delimiter='\t')

                reader.next()

                for row in reader:
                    results.append((row[0], row[1], row[2]))
                f.close()

    return render_template('index.html', results=results)


@app.route('/save', methods=['POST'])
def save():
    country = request.form['country']
    disease = request.form['disease']
    percentage = float(request.form['percentage'])

    countries = dict()
    node_label = str(uuid.uuid4())

    gm_client = Client([gearman_url])

    jobs = []
    records = 100000
    retmax = 10000

    for i in range(0, records, retmax):
        jobs.append({'task': 'process_results', 'data': {
            'disease': disease,
            'node_label': node_label,
            'retstart': i,
            'retmax': retmax}})

    completed_jobs = gm_client.submit_multiple_jobs(jobs)

    for job in completed_jobs:
        if job.result is not None:
            countries.update(job.result)

    driver = GraphDatabase.driver(
        neo4j_url,
        auth=basic_auth(neo4j_user, neo4j_password))
    session = driver.session()
    g = Graph(directed=False)

    results = session.run("MATCH (p1:`{node}`)-[:KNOWS]-(p2:`{node}`)"
                          " WITH p1, p2, count(p2) AS weight "
                          "RETURN p1.name, p2.name, weight".
                          format(node=node_label))

    for record in results:
        v1 = unescape(record['p1.name'])
        v2 = unescape(record['p2.name'])

        g.add_vertices([v1, v2])
        g.add_edge(v1, v2, weight=record['weight'])

    session.close()

    jobs = []

    for i in range(0, records, retmax):
        jobs.append({'task': 'delete_data', 'data': {
            'node_label': node_label,
            'limit': retmax}})

    gm_client.submit_multiple_jobs(jobs)

    # read network
    # extract largest connected component
    components = g.components(WEAK)
    g = components.giant()
    vs = VertexSeq(g)

    n_kols = int(round(percentage * g.vcount() / 100.0))

    # calculate centrality measures
    degrees = dict(zip(vs['name'], g.degree()))
    closeness = dict(zip(vs['name'], g.closeness()))
    eigenvector_centrality = dict(zip(vs['name'], g.eigenvector_centrality()))
    betweenness = dict(zip(vs['name'], g.betweenness()))

    # get top 15 KOLs
    top_degree = sorted(degrees.items(),
                        key=operator.itemgetter(1),
                        reverse=True)[:n_kols]
    top_closeness = sorted(closeness.items(),
                           key=operator.itemgetter(1),
                           reverse=True)[:n_kols]
    top_eigenvector = sorted(eigenvector_centrality.items(),
                             key=operator.itemgetter(1),
                             reverse=True)[:n_kols]
    top_betweenness = sorted(betweenness.items(),
                             key=operator.itemgetter(1),
                             reverse=True)[:n_kols]

    # get kols
    kols = set()
    for i in range(n_kols):
        kols.add(top_degree[i][0])
        kols.add(top_closeness[i][0])
        kols.add(top_eigenvector[i][0])
        kols.add(top_betweenness[i][0])

    # write to csv file
    us_states = {
        'AK': 'Alaska',
        'AL': 'Alabama',
        'AR': 'Arkansas',
        'AS': 'American Samoa',
        'AZ': 'Arizona',
        'CA': 'California',
        'CO': 'Colorado',
        'CT': 'Connecticut',
        'DC': 'District of Columbia',
        'DE': 'Delaware',
        'FL': 'Florida',
        'GA': 'Georgia',
        'GU': 'Guam',
        'HI': 'Hawaii',
        'IA': 'Iowa',
        'ID': 'Idaho',
        'IL': 'Illinois',
        'IN': 'Indiana',
        'KS': 'Kansas',
        'KY': 'Kentucky',
        'LA': 'Louisiana',
        'MA': 'Massachusetts',
        'MD': 'Maryland',
        'ME': 'Maine',
        'MI': 'Michigan',
        'MN': 'Minnesota',
        'MIN': 'Minneapolis',  # special case
        'MO': 'Missouri',
        'MP': 'Northern Mariana Islands',
        'MS': 'Mississippi',
        'MT': 'Montana',
        'NA': 'National',
        'NC': 'North Carolina',
        'ND': 'North Dakota',
        'NE': 'Nebraska',
        'NH': 'New Hampshire',
        'NJ': 'New Jersey',
        'NM': 'New Mexico',
        'NV': 'Nevada',
        'NY': 'New York',
        'OH': 'Ohio',
        'OK': 'Oklahoma',
        'OR': 'Oregon',
        'PA': 'Pennsylvania',
        'PHI': 'Philadelphia',  # special case
        'PR': 'Puerto Rico',
        'RI': 'Rhode Island',
        'SC': 'South Carolina',
        'SD': 'South Dakota',
        'TN': 'Tennessee',
        'TX': 'Texas',
        'UT': 'Utah',
        'VA': 'Virginia',
        'VI': 'Virgin Islands',
        'VT': 'Vermont',
        'WA': 'Washington',
        'WI': 'Wisconsin',
        'WV': 'West Virginia',
        'WY': 'Wyoming'
    }

    filename = make_safe("{}_{}_{}_{}".
                         format(
                             country,
                             disease,
                             percentage,
                             date.today().strftime('%Y-%m-%d-%H-%M-%S')))

    # print to file
    with open("{}/{}.csv".format(data_dir, filename), 'w') as f:
        f.write("FirstName\tMiddleName\tLastName\n")
        for k in kols:
            if k in countries:
                c = clean_country(countries[k])
                if c in us_states.keys() or c in us_states.values():
                    c = 'USA'

                if c.lower() == country.lower():
                    full_name = k.split('_')
                    if len(full_name) == 3:
                        f.write("{}\t{}\t{}\n".format(
                            full_name[0], full_name[1], full_name[2]))
                    else:
                        f.write("{}\t \t{}\n".format(full_name[0], full_name[1]))
        f.close()

    return redirect("/?f={}".format(filename))

if __name__ == '__main__':
    CGIHandler().run(app)