#!/usr/env/python

from __future__ import division, print_function

## Import General Tools
import sys
import os
import argparse
import logging
from datetime import datetime as dt
from datetime import timedelta as tdelta
from argparse import ArgumentParser
import re
import glob

import pymongo
from pymongo import MongoClient

from tornado.ioloop import IOLoop
from tornado.web import RequestHandler, Application, url, StaticFileHandler
from tornado import websocket
import tornado.log as tlog

from astropy import units as u
from astropy.coordinates import SkyCoord
import ephem

import IQMon

class MyStaticFileHandler(StaticFileHandler):
    def set_extra_headers(self, path):
        # Disable cache
        self.set_header('Cache-Control', 'no-store, no-cache, must-revalidate, max-age=0')


##-----------------------------------------------------------------------------
## Handler for list of images
##-----------------------------------------------------------------------------
class ListOfImages(RequestHandler):
    def get(self, telescope, subject):
        tlog.app_log.info('Get request for ListOfImages recieved')

        ## Create Telescope Object
        tlog.app_log.info('  Creating telescope object')
        config_file = os.path.join(os.path.expanduser('~'), '.{}.yaml'.format(telescope))
        tel = IQMon.Telescope(config_file)
        telescopename = tel.name
        tlog.app_log.info('  Done.')

        tlog.app_log.info('  Linking to mongo')
        client = MongoClient(tel.mongo_address, tel.mongo_port)
        tlog.app_log.info('  Connected to client.')
        db = client[tel.mongo_db]
        collection = db[tel.mongo_collection]
        tlog.app_log.info('  Retrieved collection.')

        tlog.app_log.info('  Getting list of images from mongo')

        ##---------------------------------------------------------------------
        ## If subject is formatted like a date, then get images from a date
        ##---------------------------------------------------------------------
        if re.match('\d{8}UT', subject):
            image_list = [entry for entry in\
                          collection.find({"date": subject}).sort(\
                          [('time', pymongo.ASCENDING)])]
            tlog.app_log.info('  Got list of {} images for night.'.format(len(image_list)))
        ##---------------------------------------------------------------------
        ## If subject matches a target name, then get images for that target
        ##---------------------------------------------------------------------
        else:
            tlog.app_log.info('    Getting list of target names from mongo')
            target_name_list = sorted(collection.distinct("target name"))
            if subject in target_name_list:
                tlog.app_log.info('    Getting list of image list for {} from mongo'.format(subject))
                image_list = [entry for entry in\
                              collection.find({"target name":subject}).sort(\
                              [('date', pymongo.DESCENDING),\
                              ('time', pymongo.DESCENDING)])]
                tlog.app_log.info('  Got list of {} images for target.'.format(len(image_list)))
        ##---------------------------------------------------------------------
        ## If subject is not a date or target, then render a list of targets
        ##---------------------------------------------------------------------
            else:
                image_list = []
                self.write('<html><head><style>')
                self.write('table{border-collapse:collapse;margin-left:auto;margin-right:auto;}')
                self.write('table,th,td{border:1px solid black;vertical-align:top;text-align:left;')
                self.write('padding-top:5px;padding-right:5px;padding-bottom:5px;padding-left:5px;}')
                self.write('</style></head>')
                if (len(subject) > 0) and not re.match('[tT]argets', subject):
                    self.write('<p style="text-align:center;">Could not find {} in target list:</p>'.format(subject))
                self.write('<table style="border:1px solid black;">')
                self.write('<tr><th>Target</th><th>n Images</th>')
                for target in target_name_list:
                    target_images = [entry for entry in collection.find( { "target name": target } ) ]
                    self.write('<tr><td><a href="{0}">{0}</a></td><td>{1:d}</td></tr>'.format(target, len(target_images)))
                self.write('</table></html>')

        if tel.units_for_FWHM == u.arcsec:
            FWHM_multiplier = tel.pixel_scale.value
        elif tel.units_for_FWHM == u.pix:
            FWHM_multiplier = 1.0
        else:
            FWHM_multiplier = 1.0

        if len(image_list) > 0:
            tlog.app_log.info('  Rendering ListOfImages')
            self.render("image_list.html", title="{} Results".format(telescopename),\
                        telescope = telescope,\
                        telescopename = telescopename,\
                        subject = subject,\
                        image_list = image_list,\
                        FWHM_units = tel.units_for_FWHM.to_string(),\
                        FWHM_multiplier = FWHM_multiplier,\
                       )
            tlog.app_log.info('  Done.')

##-----------------------------------------------------------------------------
## Handler for list of nights
##-----------------------------------------------------------------------------
class ListOfNights(RequestHandler):

    def get(self, telescope):
        tlog.app_log.info('Get request for ListOfNights recieved')
        telescope = telescope.strip('/')

        ## Create Telescope Object
        config_file = os.path.join(os.path.expanduser('~'), '.{}.yaml'.format(telescope))
        tel = IQMon.Telescope(config_file)
        telescopename = tel.name

        client = MongoClient(tel.mongo_address, tel.mongo_port)
        db = client[tel.mongo_db]
        collection = db[tel.mongo_collection]

        first_date_string = sorted(collection.distinct("date"), reverse=False)[0]
        first_date = dt.strptime('{} 00:00:00'.format(first_date_string), '%Y%m%dUT %H:%M:%S')
        oneday = tdelta(1, 0)
        
        tlog.app_log.info('  Building date_list')
        date_list = []
        thisdate = dt.utcnow()
        while thisdate >= first_date:
            date_list.append(thisdate.strftime('%Y%m%dUT'))
            thisdate -= oneday
        tlog.app_log.info('  Done')

        night_plot_path = os.path.abspath('/var/www/nights/')

        tlog.app_log.info('  Looping over date_list')
        nights = []
        for date_string in date_list:
            night_info = {'date': date_string }

            night_graph_file = '{}_{}.png'.format(date_string, telescope)
            if os.path.exists(os.path.join(night_plot_path, night_graph_file)):
                night_info['night graph'] = night_graph_file

            night_info['n images'] = collection.find( {"date":date_string} ).count()
            
            nights.append(night_info)
        tlog.app_log.info('  Done')

        tlog.app_log.info('  Rendering ListOfNights')
        self.render("night_list.html", title="{} Results".format(telescopename),\
                    telescope = telescope,\
                    telescopename = telescopename,\
                    nights = nights,\
                   )
        tlog.app_log.info('  Done')


##-----------------------------------------------------------------------------
## Main
##-----------------------------------------------------------------------------
def main():
    ##-------------------------------------------------------------------------
    ## Parse Command Line Arguments
    ##-------------------------------------------------------------------------
    ## create a parser object for understanding command-line arguments
    parser = ArgumentParser(description="Describe the script")
    ## add flags
    parser.add_argument("--with-status",
        action="store_true", dest="status",
        default=False, help="Use status handler from custom_handlers.py")
    args = parser.parse_args()

    LogConsoleHandler = logging.StreamHandler()
    LogConsoleHandler.setLevel(logging.DEBUG)
    LogFormat = logging.Formatter('%(asctime)23s %(levelname)8s: %(message)s')
    LogConsoleHandler.setFormatter(LogFormat)
    tlog.app_log.addHandler(LogConsoleHandler)
    tlog.app_log.setLevel(logging.DEBUG)
    tlog.access_log.addHandler(LogConsoleHandler)
    tlog.access_log.setLevel(logging.DEBUG)
    tlog.gen_log.addHandler(LogConsoleHandler)
    tlog.gen_log.setLevel(logging.DEBUG)

    list_of_handlers = [
                        url(r"/(\w+/?$)", ListOfNights),
                        url(r"/(\w+)/(\w+)", ListOfImages),
                        (r"/static/(.*)", MyStaticFileHandler, {"path": "/var/www"}),
                       ]

    if args.status:
        tlog.app_log.info('Importing status handler')
        try:
            from custom_handlers import Status
            list_of_handlers.append(url(r"/", Status))
            tlog.app_log.info('  Done')
        except:
            tlog.app_log.warning('  Failed')
            pass


    app = Application(list_of_handlers)
    app.listen(80)
    IOLoop.current().start()

if __name__ == '__main__':
    main()
