#! /usr/bin/env python

import datetime
import logging
import time

from Acspy.Clients.SimpleClient import PySimpleClient
from Acspy.Util.ACSCorba import getManager

logging.basicConfig(
    filename='/home/gavino/logs/receivers.log',
    format='%(asctime)s\t%(message)s')

common_receiver_properties = (
    'cryoTemperatureCoolHead', 
    'cryoTemperatureCoolHeadWindow', 
    'cryoTemperatureLNA', 
    'cryoTemperatureLNAWindow', 
    'environmentTemperature',
    'vacuum')

mapping = (
    ('ANTENNA/Boss', ('rawAzimuth', 'rawElevation')),
    ('RECEIVERS/SRTLPBandReceiver', common_receiver_properties),
    ('RECEIVERS/SRTKBandMFReceiver', common_receiver_properties),
    ('RECEIVERS/SRT7GHzReceiver', common_receiver_properties),
)


timing = 180 #seconds

while True:
    if getManager():
        try:
            components = []
            for component_name, property_names in mapping:
                client = PySimpleClient()
                try:
                    component = client.getComponent(component_name)
                    components.append((component, property_names))
                except:
                    pass

            for component, property_names in components:
                for pname in property_names:
                    get_property_obj = getattr(component, '_get_%s' % pname)
                    property_obj = get_property_obj()
                    value, completion = property_obj.get_sync()
                    cname = component._get_name()
                    t = datetime.datetime.now()
                    line = '%s.%s' % (cname, pname)
                    logging.info('  %s%e' % (line.ljust(65), value))

        except KeyboardInterrupt:
            logging.info('program closed by the user')
            raise
        except Exception, ex:
            pass
        finally:
            try:
                if getManager():
                    client.disconnect()
            except:
                logging.error('can not disconnect the client')

    time.sleep(timing)
