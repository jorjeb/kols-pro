# -*- coding: utf-8 -*-
from hashids import Hashids
from flask import Flask, render_template, request, redirect
from redis import Redis
from rq import Queue, cancel_job
import csv
import json
import kols
import os
import time

app = Flask(__name__)
q = Queue(connection=Redis())

data_dir = 'data'  # no trailing slash


@app.route('/', defaults={'f': None})
@app.route('/<string:f>')
def index(f):
    run_checker = False
    job_id = request.args['job'] if 'job' in request.args else False
    results = []

    if 'check' in request.args:
        if f is not None:
            if os.path.isfile("{}/.{}".format(data_dir, f)):
                return json.dumps({'done': False})

        return json.dumps({'done': True})

    if f is not None:
        if os.path.isfile("{}/.{}".format(data_dir, f)):
            run_checker = True

        filename = "{}/{}.csv".format(data_dir, f)

        if os.path.isfile(filename):
            with open(filename, 'r') as f:
                reader = csv.reader(f, delimiter='\t')

                next(reader)

                for row in reader:
                    results.append((row[0], row[1], row[2]))
                f.close()

    return render_template('index.html',
                           results=results,
                           f=f,
                           run_checker=run_checker,
                           job_id=job_id)


@app.route('/cancel')
def cancel():
    if 'job' in request.args:
        cancel_job(request.args['job'], connection=Redis())

    return redirect('/kols')


@app.route('/save', methods=['POST'])
def save():
    country = request.form['country']
    disease = request.form['disease']
    percentage = float(request.form['percentage'])
    filename = Hashids().encode(int(time.time()))

    open("{}/.{}".format(data_dir, filename), 'a').close()

    job = q.enqueue(kols.save, country, disease, percentage, filename, timeout=1800)

    return redirect("/kols/{}?job={}".format(filename, job.id))

if __name__ == '__main__':
    app.run(debug=True)
