import pymongo

from iqmon import get_webpage_config


def insert_mongodoc(collection, mongodoc, log=None):
    cfg = get_webpage_config()
    if log: log.info(f'Connecting to mongoDB')
    mongo_host = cfg['mongo'].get('host')
    mongo_port = cfg['mongo'].getint('port')
    mongo_db = cfg['mongo'].get('db')
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


def query_non_weather_devices(cfg=None):
    if cfg is None:
        cfg = get_config()

    # Telescope
    if cfg['devices'].get('telescope', None) is not None:
        query_type, 
