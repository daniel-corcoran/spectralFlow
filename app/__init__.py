from flask import Flask


template_session = {'samples': [],
                    'query': [],
                    'date': None,
                    'cosine_threshold': 0,
                    'abundance_threshold': 0}
active_sessions = {}
# key = token
# value = files, settings, etc.


app = Flask(__name__)

