import sys
import pymongo
import time

from iqmon import get_webpage_config

def insert_mongodoc(collection, mongodoc, log=None):
    cfg = get_webpage_config()
    mongo_host = cfg['mongo'].get('host')
    mongo_port = cfg['mongo'].getint('port')
    mongo_db = cfg['mongo'].get('db')
    if log: log.info(f'Connecting to mongo: {mongo_host}:{mongo_port} {mongo_db}.{collection}')
    mongoclient = pymongo.MongoClient(mongo_host, mongo_port)
    collection = mongoclient[mongo_db][collection]
    try:
        inserted_id = collection.insert_one(mongodoc).inserted_id
        log.info("  Inserted document with id: {}".format(inserted_id))
    except:
        e = sys.exc_info()[0]
        if log: log.error('Failed to add new document')
        if log: log.error(e)
        else: print(e)
    mongoclient.close()
