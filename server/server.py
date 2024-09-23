#! /usr/bin/env python

import os
import uuid
import re
import subprocess
import argparse
import json
from subprocess import call
from flask import Flask, send_file, flash, send_from_directory, request, redirect, url_for, jsonify
from flask_cors import CORS
from werkzeug.utils import secure_filename

WALLYWS = os.path.dirname(os.path.abspath(__file__))

app = Flask(__name__)
CORS(app)
app.config['WALLY'] = os.path.join(WALLYWS, "..")
app.config['UPLOAD_FOLDER'] = os.path.join(app.config['WALLY'], "data")
app.config['MAX_CONTENT_LENGTH'] = 8 * 1024 * 1024   #maximum of 8MB

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in set(['png'])

uuid_re = re.compile(r'^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}-zoom[0-9]*$')

def is_valid_uuid(s):
    return uuid_re.match(s) is not None

@app.route('/api/v1/download/<uuid>/<ext>')
def download(uuid, ext):
    if is_valid_uuid(uuid):
        filename = "" + uuid + "." + ext
        if allowed_file(filename):
            sf = os.path.join(app.config['UPLOAD_FOLDER'], uuid[0:2])
            if os.path.exists(sf):
                if os.path.isfile(os.path.join(sf, filename)):
                    return send_file(os.path.join(sf, filename), mimetype='image/gif')
    return "File does not exist!"

@app.route('/api/v1/upload', methods=['POST'])
def upload_file():
    if request.method == 'POST':
        uuidstr = str(uuid.uuid4())
        
        # Get subfolder
        sf = os.path.join(app.config['UPLOAD_FOLDER'], uuidstr[0:2])
        if not os.path.exists(sf):
            os.makedirs(sf)

        ds = None
        if 'dataset' in request.form.keys():
            ds = request.form['dataset']
        else:
            return jsonify(errors = [{"title": "Data set is missing!"}]), 400
        if ds == "1000 Genomes ONT Vienna GRCh38/hg38":
            ds = "1kGP_ONT_hg38"
        elif ds == "1000 Genomes ONT Vienna T2T-CHM13":
            ds = "1kGP_ONT_t2t"
        elif ds == "1000 Genomes Illumina high-coverage GRCh38/hg38":
            ds = "1kGP_ILL_hg38"
        else:
            return jsonify(errors = [{"title": "Data set does not exist!"}]), 400
        chrName = None
        if 'chr' in request.form.keys():
            chrName = request.form['chr']
        else:
            return jsonify(errors = [{"title": "Chromosome is missing!"}]), 400
        regionStart = 0
        if 'regionStart' in request.form.keys():
            regionStart = int(request.form['regionStart'])
        else:
            return jsonify(errors = [{"title": "Region start is incorrect!"}]), 400
        regionEnd = 0
        if 'regionEnd' in request.form.keys():
            regionEnd = int(request.form['regionEnd'])
        else:
            return jsonify(errors = [{"title": "Region end is incorrect!"}]), 400
        if regionStart >= regionEnd:
            return jsonify(errors = [{"title": "Region start needs to be smaller than region end!"}]), 400
        if (regionEnd - regionStart) >= 100000:
            return jsonify(errors = [{"title": "Region needs to be smaller than 100Kbp!"}]), 400
        if (regionEnd - regionStart) <= 10:
            return jsonify(errors = [{"title": "Region needs to be larger than 10bp!"}]), 400
        samples = None
        if 'samples' in request.form.keys():
            samples = request.form['samples']
            s = set(samples.replace('\r\n','\n').strip().split('\n'))
            if len(s) < 1:
                return jsonify(errors = [{"title": "Please provide sample names!"}]), 400
            if len(s) > 20:
                return jsonify(errors = [{"title": "Maximum number of samples is 20!"}]), 400
            samples = ' '.join(s)
            if samples == '':
                return jsonify(errors = [{"title": "Please provide sample names!"}]), 400
        
        # Run wally
        outdir = os.path.join(sf)
        logfile = os.path.join(sf, "wally_" + uuidstr + ".log")
        errfile = os.path.join(sf, "wally_" + uuidstr + ".err")
        blexe = os.path.join(app.config['WALLY'], "wallyapp.sh")
        reg = chrName + ":" + str(regionStart) + "-" + str(regionEnd) + ":" + uuidstr
        with open(logfile, "w") as log:
            with open(errfile, "w") as err:
                try:
                    return_code = call([blexe, outdir, ds, reg, samples], stdout=log, stderr=err)
                except OSError as e:
                    if e.errno == os.errno.ENOENT:
                        return jsonify(errors = [{"title": "WallyApp script not found!"}]), 400
                    else:
                        return jsonify(errors = [{"title": "OSError " + str(e.errno)  + " running WallyApp script!"}]), 400
        if return_code != 0:
            errInfo = "!"
            with open(errfile, "r") as err:
                errInfo = ": " + err.read()
            return jsonify(errors = [{"title": "Error in running Wally" + errInfo}]), 400
        urlout = "download/" + uuidstr + "-zoom5"
        dt = {}
        dt["url"] = urlout
        return jsonify(data=dt), 200
    return jsonify(errors = [{"title": "Error in handling POST request!"}]), 400

@app.route('/api/v1/health', methods=['GET'])
def health():
    return jsonify(status="OK")

if __name__ == '__main__':
    app.run(host = '0.0.0.0', port=3300, debug = True, threaded=True)
