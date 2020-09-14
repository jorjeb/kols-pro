# -*- coding: utf-8 -*-
from Bio import Entrez
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor
from igraph import Graph, operator, VertexSeq, WEAK
import os
import time
import traceback

data_dir = 'data'  # no trailing slash


def search_field(disease, retstart, retmax):
    Entrez.email = ''
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retstart=retstart,
                            retmax=retmax,
                            retmode='xml',
                            term=disease)
    results = Entrez.read(handle)
    return results


def fetch_details(id_list):
    if not id_list:
        return dict()

    ids = ','.join(id_list)
    Entrez.email = ''
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results


def find_between_r(s, first, last):
    try:
        start = s.rindex(first) + len(first) + 1
        end = s.rindex(last, start)
        return s[start:end]
    except ValueError:
        return None


def process_results(data):
    try:
        disease = data['disease']
        retstart = data['retstart']
        retmax = data['retmax']

        print("Start fetching records {} to {} ...".format(retstart + 1,
                                                           retstart + retmax))

        results = search_field("\"{}\"".format(disease), retstart, retmax)
        id_list = results['IdList']

        papers = fetch_details(id_list)
        countries = dict()
        edges = defaultdict(list)

        for i, paper in enumerate(papers):
            if 'MedlineCitation' in paper.keys():
                if 'Article' in paper['MedlineCitation'].keys():
                    if 'AuthorList' in paper['MedlineCitation']['Article'].keys():
                        author_list = (paper['MedlineCitation']['Article']
                                       ['AuthorList'])
                        authors = set()
                        for index in range(0, min(len(author_list), 12)):
                            a = author_list[index]
                            name = ''
                            if 'LastName' in a.keys() and 'ForeName' in a.keys():
                                name = "{}_{}".format(
                                    a['ForeName'].replace(' ', '_'),
                                    a['LastName'].replace(' ', ''))

                                authors.add(name)
                            if (name not in countries.keys()):
                                if 'AffiliationInfo' in a.keys():
                                    if len(a['AffiliationInfo']) > 0:
                                        if ('Affiliation' in
                                                a['AffiliationInfo'][0].keys()):
                                            s_country = (find_between_r(
                                                (a['AffiliationInfo'][0]
                                                    ['Affiliation']),
                                                ',', '.'))
                                            if s_country is not None:
                                                if 'China' in s_country:
                                                    s_country = 'China'
                                                if 'United States' in s_country:
                                                    s_country = 'USA'
                                                countries[name] = s_country

                        authors = list(authors)

                        # create edges
                        for i in range(len(authors)):
                            u = authors[i]

                            for j in range(i + 1, len(authors)):
                                v = authors[j]

                                edges[u].append(v)
                                edges[v].append(u)

        print("Done fetching records {} to {} ...".format(retstart + 1, retstart + retmax))

        return (edges, countries)
    except:
        traceback.print_exc()


def clean_country(c):
    try:
        pos = c.rindex('(')
        return c[0:pos - 1]
    except ValueError:
        return c


def save(country, disease, percentage, filename):
    percentage = float(percentage)

    countries = dict()
    records = 1000
    retmax = 100

    tuple_list = []

    with ProcessPoolExecutor() as pool:
        jobs = []
        jobs_done = []

        for i in range(0, records, retmax):
            job = pool.submit(process_results, {
                'disease': disease,
                'retstart': i,
                'retmax': retmax})

            setattr(job, 'tag', i)

            jobs.append(job)

        while len(jobs_done) < len(jobs):
            for job in jobs:
                tag = getattr(job, 'tag')

                if job.done() and tag not in jobs_done:
                    jobs_done.append(tag)

            time.sleep(1)

        network_edges = defaultdict(list)

        for job in jobs:
            result = job.result()

            if result is not None:
                edges = result[0]
                info = result[1]

                if info is not None:
                    countries.update(info)

                if edges is not None:
                    for k in edges.keys():
                        if k in network_edges:
                            network_edges[k] += edges[k]
                        else:
                            network_edges[k] = edges[k]

        for k in network_edges.keys():
            neighbors = Counter(network_edges[k])

            for n in neighbors.keys():
                if n > k:
                    tuple_list.append((n, k, int(neighbors[n])))

    # read network
    # extract largest connected component
    g = Graph.TupleList(tuple_list, weights=True)
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

    os.remove("{}/.{}".format(data_dir, filename))
