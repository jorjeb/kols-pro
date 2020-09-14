from Bio import Entrez
from neo4j.v1 import GraphDatabase, basic_auth
from py2neo.database.cypher import cypher_escape
from sets import Set
import gearman
import json
import traceback

gearman_url = '127.0.0.1:4730'
neo4j_url = 'bolt://localhost'
neo4j_user = 'neo4j'
neo4j_password = '1234'


class DataEncoder(gearman.DataEncoder):
    @classmethod
    def encode(cls, object):
        return json.dumps(object)

    @classmethod
    def decode(cls, json_string):
        return json.loads(json_string)


class Worker(gearman.GearmanWorker):
    data_encoder = DataEncoder

    def on_job_exception(self, current_job, exc_info):
        print(exc_info)
        print(''.join(traceback.format_tb(exc_info[2])))

        return super(Worker, self).on_job_exception(current_job, exc_info)


def search_field(disease, restart, retmax):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retstart=restart,
                            retmax=retmax,
                            retmode='xml',
                            term=disease)
    results = Entrez.read(handle)
    return results


def fetch_details(id_list):
    if not id_list:
        return dict()

    ids = ','.join(id_list)
    Entrez.email = 'your.email@example.com'
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
        return "Not found"


def delete_data(gearman_worker, gearman_job):
    node_label = gearman_job.data['node_label']
    limit = gearman_job.data['limit']

    driver = GraphDatabase.driver(
        neo4j_url,
        auth=basic_auth(neo4j_user, neo4j_password))
    session = driver.session()
    session.run("MATCH (n:`{}`) WITH n LIMIT {} DETACH DELETE n".
                format(node_label, limit))

    session.close()


def process_results(gearman_worker, gearman_job):
    disease = gearman_job.data['disease']
    node_label = gearman_job.data['node_label']
    retstart = gearman_job.data['retstart']
    retmax = gearman_job.data['retmax']

    results = search_field("\"{}\"".format(disease), retstart, retmax)
    id_list = results['IdList']
    papers = fetch_details(id_list)

    countries = dict()

    driver = GraphDatabase.driver(
        neo4j_url,
        auth=basic_auth(neo4j_user, neo4j_password))
    session = driver.session()

    for i, paper in enumerate(papers):
        if 'MedlineCitation' in paper.keys():
            if 'Article' in paper['MedlineCitation'].keys():
                if 'AuthorList' in paper['MedlineCitation']['Article'].keys():
                    author_list = (paper['MedlineCitation']['Article']
                                   ['AuthorList'])
                    authors = Set()
                    for index in range(0, min(len(author_list), 12)):
                        a = author_list[index]
                        name = ''
                        if 'LastName' in a.keys() and 'ForeName' in a.keys():
                            name = "{}_{}".format(
                                a['ForeName'].encode('utf-8').replace(' ', '_'),
                                a['LastName'].encode('utf-8').replace(' ', ''))

                            authors.add(name)
                        if (name not in countries.keys() or
                                countries[name] == "Not found"):
                            if 'AffiliationInfo' in a.keys():
                                if len(a['AffiliationInfo']) > 0:
                                    if ('Affiliation' in
                                            a['AffiliationInfo'][0].keys()):
                                        s_country = (find_between_r(
                                            (a['AffiliationInfo'][0]
                                                ['Affiliation']),
                                            ',', '.'))
                                        if 'China' in s_country:
                                            s_country = 'China'
                                        if 'United States' in s_country:
                                            s_country = 'USA'
                                        countries[name] = s_country
                                else:
                                    countries[name] = "Not found"

                    authors = list(authors)

                    # create network
                    for i in range(len(authors)):
                        u = authors[i]

                        for j in range(i + 1, len(authors)):
                            v = authors[j]
                            session.run("MERGE (p1:`{node}` {{name:\"{}\"}}) "
                                        "MERGE (p2:`{node}` {{name:\"{}\"}}) "
                                        "CREATE (p1)-[:KNOWS]->(p2)".
                                        format(cypher_escape(u).encode('utf-8'),
                                               cypher_escape(v).encode('utf-8'),
                                               node=node_label))

    session.close()

    return countries

if __name__ == '__main__':
    gm_worker = Worker([gearman_url])
    gm_worker.register_task('process_results', process_results)
    gm_worker.register_task('delete_data', delete_data)
    gm_worker.work()
