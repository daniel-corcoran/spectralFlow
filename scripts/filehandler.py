from app import app
from flask import request
from flask import render_template
import random
import string
import os
import json

UPLOAD_FOLDER = 'files/'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
ALLOWED_EXTENSIONS = {'mgf', 'MGF'}


def get_random_string(length):
    letters_and_digits = string.ascii_letters + string.digits
    return ''.join((random.choice(letters_and_digits) for i in range(length)))


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/API/upload', methods=['POST', 'GET'])
def upload():

    # Error types:
    """
    Error types:
    err_1: No reference file uploaded
    err_2: No sample file uploaded

    err_3: reference file bad extension or filename
    err_4: sample file bad extension or filename

    err_5: reference file exceeds file size limit (25mb?)
    err_6: sample file exceeds file size limit (25mb?)

    err_7: reference file could not be interpreted
    err_8: sample file could not be interpreted
    """

    response_json = {'success': True,
                     'session_token': 0,
                     'err_1': False,
                     'err_2': False,
                     'err_3': False,
                     'err_4': False,
                     'err_5': False,
                     'err_6': False,
                     'err_7': False,
                     'err_8': False}

    if request.method == 'POST':

        if 'reference' not in request.files:

            response_json['success'] = False
            response_json['err_1'] = True

        if 'sample' not in request.files:

            response_json['success'] = False
            response_json['err_2'] = True

        if 'reference' in request.files and 'sample' in request.files:

            reference = request.files['reference']
            sample = request.files['sample']

            if reference.filename == '':

                response_json['success'] = False
                response_json['err_1'] = True

            if sample.filename == '':

                response_json['success'] = False
                response_json['err_2'] = True

            if not response_json['success']:
                # Stop processing here because we can't interpret the files.
                # Return error message back to request.

                return json.dumps(response_json)

            else:
                # So far so good, we can continue.
                tkn = get_random_string(128)

                if reference and allowed_file(reference.filename):

                    filename = tkn + 'reference.mgf'
                    reference.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

                else:

                    response_json['success'] = False
                    response_json['err_3'] = True

                if sample and allowed_file(sample.filename):

                    filename = tkn + 'sample.mgf'
                    sample.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

                else:
                    response_json['success'] = False
                    response_json['err_4'] = True

                if response_json['success']:

                    # If we've reached this point, it means that our files at least both exist and have the correct extension.
                    # Now we need to see if they meet the size and formatting requirements.
                    # FIXME: Check to see if files are readable and below size restrictions.

                    response_json['session_token'] = tkn
                    return json.dumps(response_json)

        return json.dumps(response_json)

    else:
        return render_template('start.html')


