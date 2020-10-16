from flask import Flask
from flask import render_template
from app import app

from scripts import filehandler, routes, interfaceAPI

import waitress



if __name__ == '__main__':
    try:
        waitress.serve(app, host='0.0.0.0', port=80)
    except:
        waitress.serve(app, host='0.0.0.0', port=8000)
