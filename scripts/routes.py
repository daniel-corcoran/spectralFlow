# General routing features

from app import app
from flask import render_template, request





@app.route('/interface')
def interface():

    # Once the user has sucessfully uploaded the files, we load the main interface where we perform spectral search.

    token = request.args.get('session_id', type=str, default=-1)

    print(token)

    return render_template('interface.html', token = token)




@app.route('/')
def entrypoint():
    return render_template('start.html')

