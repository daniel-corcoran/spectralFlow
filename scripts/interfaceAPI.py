from app import app
from flask import request, render_template, Request
from search.spectralFlow.flow import yield_vector_from_file
import json
import matplotlib

import matplotlib.pyplot as plt
import io
import base64
import time
import random
import os
from search.spectralFlow.flow import web_api_hook

avail = True


def render_table(response_json, mode):
        num = len(response_json)
        print("I have to process {} spectra.".format(num))
        stime = time.time()
        for a, b in enumerate(response_json):
            response_json[a]['pepmass'] = round(response_json[a]['pepmass'], 4)


        return render_template('iontable.html', data=[[a, b] for a, b in enumerate(response_json)], mode=mode)


@app.route('/API/search', methods=['POST', 'GET'])
def request_spectral_match():
    # TODO: Can we show results as they load in?
    mz = {'on': True, None: False}[request.form.get('mz')]
    cos = request.form.get('cos', default=0, type=float)
    ab = request.form.get('ab', default=0.07, type=float)
    tkn = request.form.get('tkn', default=-1, type=str)
    print(mz)
    print(cos)
    print(ab)
    print(tkn)

    response_list = web_api_hook(ses_tkn = tkn,
                 cos_threshold = cos,
                 mass_check = mz,
                 mass_delta = 0.004, # FIXME
                 abundance_cut=ab)
    # Start a thread for spectral matching.
    results_table = render_template('resultstable.html', results=response_list)
    # Construct a new HTML table with the results.


    return {'success': True, 'html': results_table}

@app.route("/API/mass_spectrum", methods=['POST'])
def render_mass_spectrum():

    try:
        content = request.json
        mz = content['mass']
        ab = content['abundance']

        fig = plt.figure(figsize = (5,2))
        ax = fig.add_subplot(1,1,1)
        ax.set_yticks([t / 5 for t in range(6)])
        ax.set_xticks([round(min(mz), 0), round(max(mz), 0)])
        ax.set_ylim(0, 1.2)
        for a, b in zip(mz, ab):
            ax.axvline(x=a, ymin=0, ymax=b * 0.83)


        fid = random.randint(1, 10000)
        fig.savefig(f"{fid}.png", transparent=False, dpi=50)
        with open(f"{fid}.png", 'rb') as plot_file:
            base64_string = base64.b64encode(plot_file.read()).decode()
            img_elem = '<img src="data:image/png;base64,{}"\>'.format(base64_string)

        plt.close(fig)
        os.system('rm ' + f"{fid}.png")
        return json.dumps({'img': img_elem})
    except Exception as E:
        print(E)
        json.dumps({'img': "Could not load spectrum"})


@app.route('/API/reference_ion_json')
def reference_ion_json():
    # Returns the list of reference ions for a session key
    session_id = request.args.get('session_id', type=str, default='')
    ab_cut = request.args.get('relative_ion_cutoff', type=float, default=0)

    gen = yield_vector_from_file(f'files/{session_id}reference.mgf', ab_cut=ab_cut)
    gen = [x for x in gen]
    response_json = {}

    for count, data in enumerate(gen):
        response_json[count] = data
    response_json = {'ions': response_json, 'html': render_table(gen, 'reference')}
    return json.dumps(response_json)


@app.route('/API/sample_ion_json')
def sample_ion_json():
    # Returns the list of reference ions for a session key
    session_id = request.args.get('session_id', type=str, default='')
    ab_cut = request.args.get('relative_ion_cutoff', type=float, default=0)

    gen = yield_vector_from_file(f'files/{session_id}sample.mgf', ab_cut=ab_cut)
    gen = [x for x in gen]

    response_json = {}

    for count, data in enumerate(gen):
        response_json[count] = data
    response_json = {'ions': response_json, 'html': render_table(gen, 'sample')}

    return json.dumps(response_json)
