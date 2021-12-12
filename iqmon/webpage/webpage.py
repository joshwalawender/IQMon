from pathlib import Path
import configparser
import logging
import flask
import pymongo
from datetime import datetime, timedelta


##-------------------------------------------------------------------------
## Create logger object
##-------------------------------------------------------------------------
log = logging.getLogger('FlaskLogger')
log.setLevel(logging.DEBUG)
LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s',
                              datefmt='%Y-%m-%d %H:%M:%S')
## Set up console output
## Set up file output
LogFileName = '/usr/local/var/log/flask.log'
LogFileHandler = logging.FileHandler(LogFileName)
LogFileHandler.setLevel(logging.DEBUG)
LogFileHandler.setFormatter(LogFormat)
log.addHandler(LogFileHandler)


app = flask.Flask(__name__)

@app.route("/")
def status():
    log.info(f'Building {__name__} status')
    content = "System status"
    return content


@app.route("/<string:telescope>/targets")
def targetList(telescope):
    log.info(f'Building {__name__} targetList {telescope}')
    template_path = Path(__file__).parent / 'targetList.html'
    return flask.render_template(template_path, telescope=telescope)


@app.route("/<string:telescope>/images/<string:date>")
def imageList(telescope, date):
    log.info(f'Building {__name__} imageList {telescope} {date}')
    subject = date

    if True:
        cfg_path = Path(__file__).parent.parent / 'configs' / 'pipeline.cfg'
        log.debug(f'Reading config file {cfg_path}')
        cfg = configparser.ConfigParser()
        cfg.read(cfg_path)
        log.debug('Connecting to mongo db')
        mongo_host = cfg['mongo'].get('host')
        log.debug(f"  mongo_host : {mongo_host}")
        mongo_port = cfg['mongo'].getint('port')
        log.debug(f"  mongo_port : {mongo_port}")
        mongo_db = cfg['mongo'].get('db')
        log.debug(f"  mongo_db : {mongo_db}")
        mongo_collection = cfg['mongo'].get('collection')
        log.debug(f"  mongo_collection : {mongo_collection}")
        mongoclient = pymongo.MongoClient(mongo_host, mongo_port)
        mongo_iqmon = mongoclient[mongo_db][mongo_collection]

        start = datetime.strptime(date, '%Y%m%dUT')
        end = start + timedelta(days=1)
        query_dict = {'telescope': telescope,
                      'date': {'$gt': start, '$lt': end}}
        image_list = [d for d in mongo_iqmon.find(query_dict,
                      sort=[('date', pymongo.ASCENDING)] ) ]
    else:
#     except Exception as e:
        log.error('Could not connect to mongo db')
        log.error(e)
        image_list = []

    log.info(f"Rendering template")
    return flask.render_template('imageList.html',
                                 telescope=telescope,
                                 subject=subject,
                                 image_list=image_list,
                                 )


if __name__ == "__main__":
    app.run()
