#!/usr/bin/env python
# Author: Marco Buttu <m.buttu@oa-cagliari.inaf.it>
# Copyright: This module has been placed in the public domain.

"""This module defines several functions that make the command answers."""

from parameters import closers, time_stamp, number_of_axis, app_nr, db_name, response_types, app_state_max, cab_state_max 
import random
import posutils

# The actual value of SRP position
srp_act_position = [1570.20]*number_of_axis['SRP']

# Drive cabinet states
dc_startup_state = 1
dc_park_state = 3
dc_ok = 0

# Application states
app_remote_auto = 4

def getpos(cmd_num, app_num, response_type):
    """This function make an answer for a getpos request.

    Parameters:

    - `cmd_num`: the command identification number
    - `app_num`: the application number (minor servo address)
    - `response_type`: the type of response: expected_ack, mixed

    Return a list of answers.
    """
    servo_type = app_nr[app_num]
    # The default response value 13477... is the time
    mixed = []
    mixed.append('?' + 'getpos' + ':0=%d> %s' % (random.randrange(0,5), time_stamp))
    mixed.append('?' + 'getpos' + ':0=%d> %s' % (app_num, time_stamp))
    mixed.append('?' + 'getpos' + ':0=#none> 1000')
    db = posutils.PositionDB()
    data = db.get(servo_type);
    expected_ack = '?' + 'getpos' + ':%d=%d> %s' % (cmd_num, app_num, posutils.now())
    # Read the positions stored in a shelve db by a setpos command
    for item in data[1:]:
        expected_ack += ',%s' %item

    return response_list('getpos', cmd_num, app_num, expected_ack, mixed, response_type)


def getappstatus(cmd_num, app_num, response_type):
    """This function make an answer for a getappstatus request.

    Parameters:

    - `cmd_num`: the command identification number
    - `app_num`: the application number (minor servo address)
    - `response_type`: the type of response: expected_ack, mixed

    Return a list of answers.
    """
    # The default response value 13477... is the time
    value = '0000030D'
    expected_ack = '?' + 'getappstatus' + ':%d=%d> %s' %(cmd_num, app_num, value)
    mixed = []
    mixed.append('?' + 'getappstatus' + ':0=%d> %s' % (random.randrange(0,5), value))
    mixed.append('?' + 'getappstatus' + ':0=%d> %s' % (app_num, value))
    mixed.append('?' + 'getappstatus' + ':0=#none> 1000')

    return response_list('getappstatus', cmd_num, app_num, expected_ack, mixed, response_type)


def getstatus(cmd_num, app_num, response_type):
    """This function make an answer for a getstatus request.

    Parameters:

    - `cmd_num`: the command identification number
    - `app_num`: the application number (minor servo address)
    - `response_type`: the type of response: expected_ack, mixed

    Return a list of answers.
    """
    servo_type = app_nr[app_num]

    mixed = []
    mixed.append('?' + 'getstatus' + ':0=%d> ' % (random.randrange(0,5),))
    mixed.append('?' + 'getstatus' + ':0=%d> ' % (app_num,))
    mixed.append('?' + 'getstatus' + ':0=#none> 1000')

    # Random state
    # app_state = random.randrange(0, app_state_max)
    # app_status = ''.join([str(item) for item in [hex(random.randrange(16)).split('x')[-1].upper() for i in range(4)]])
    # cab_state = random.randrange(0, cab_state_max)

    # Setup state
    # cab_state = dc_startup_state

    # Park state
    # cab_state = dc_park_state

    # Ready
    app_state = app_remote_auto
    app_status = "FFFF"
    cab_state = dc_ok

    expected_ack = '?' + 'getstatus' + ':%d=%d> ' % (cmd_num, app_num)
    db = posutils.PositionDB()
    data = db.get(servo_type)
    expected_ack += '%d,%d,%s,%d' %( posutils.now(), app_state, app_status, cab_state)
    # Read the positions stored in a shelve db by a setpos command
    for item in data[1:]:
        expected_ack += ',%s' %item

    return response_list('getstatus', cmd_num, app_num, expected_ack, mixed, response_type)


def setpos(cmd_num, app_num, response_type, *params):
    """This function make an answer for a setpos request.

    Parameters:

    - `cmd_num`: the command identification number
    - `app_num`: the application number (minor servo address)
    - `response_type`: the type of response: expected_ack, mixed

    Return a list of answers.
    """
    servo_type = app_nr[app_num]
    db = posutils.PositionDB()
    db.insert(servo_type, params[-number_of_axis[servo_type]:], params[0])
    expected_ack = '@' + 'setpos' + ':%d=%d' % (cmd_num, app_num)
    for param in params:
        expected_ack += ",%s" %param
    mixed =  []
    mixed.append('?' + 'setpos' + ':0=%d> 0' %random.randrange(0,5))
    mixed.append('?' + 'setpos' + ':0=%d> 0' %app_num)
    mixed.append('@' + 'setpos' + ':0=foo>')

    return response_list('setpos', cmd_num, app_num, expected_ack, mixed, response_type)


def setup(cmd_num, app_num, response_type, *params):
    """This function make an answer for a setup request.

    Parameters:

    - `cmd_num`: the command identification number
    - `app_num`: the application number (minor servo address)
    - `response_type`: the type of response: expected_ack, mixed

    Return a list of answers.
    """
    expected_ack = '@' + 'setup' + ':%d=%d' % (cmd_num, app_num)
    for param in params:
        expected_ack += ",%s" %param
    mixed = []
    mixed.append('?' + 'stow' + ':0=%d> 0' %random.randrange(0,5))
    mixed.append('?' + 'setup' + ':0=%d> 0' %app_num)

    return response_list('setup', cmd_num, app_num, expected_ack, mixed, response_type)


def stow(cmd_num, app_num, response_type, *params):
    """This function make an answer for a stow request.

    Parameters:

    - `cmd_num`: the command identification number
    - `app_num`: the application number (minor servo address)
    - `response_type`: the type of response: expected_ack, mixed

    Return a list of answers.
    """
    expected_ack = '@' + 'stow' + ':%d=%d' % (cmd_num, app_num)
    for param in params:
        expected_ack += ",%s" %param
    mixed = []
    mixed.append('?' + 'setup' + ':0=%d> 0' %random.randrange(0,5))
    mixed.append('?' + 'stow' + ':0=%d> 0' %app_num)

    return response_list('stow', cmd_num, app_num, expected_ack, mixed, response_type)


def clean(cmd_num, app_num, response_type, *params):
    """This function make an answer for a setup request.

    Parameters:

    - `cmd_num`: the command identification number
    - `app_num`: the application number (minor servo address)
    - `response_type`: the type of response: expected_ack, mixed

    Return a list of answers.
    """
    expected_ack = '@' + 'clean' + ':%d=%d' % (cmd_num, app_num)
    for param in params:
        expected_ack += ",%s" %param
    mixed = []
    mixed.append('?' + 'stow' + ':0=%d> 0' %random.randrange(0,5))
    mixed.append('?' + 'clean' + ':0=%d> 0' %app_num)

    return response_list('clean', cmd_num, app_num, expected_ack, mixed, response_type)




def getspar(cmd_num, app_num, response_type, *params):
    """This function make an answer all the getspar requestes.

    Parameters:

    - `cmd_num`: the command identification number
    - `app_num`: the application number (minor servo address)
    - `response_type`: the type of response: expected_ack, mixed
    - `params`: a tuple containing the index and sub-index of a parameter request.

    Return a list of answers.
    """
    expected_ack = '?' + 'getspar' + ':%d=%d' % (cmd_num, app_num)
    nfm = '?' + 'getspar' + ':0=%d' %app_num

    for param in params[1:]:
        expected_ack += ',%s' %param
        nfm += ',%s' %param

    expected_ack += '> %d' % (sum([ int(param) for param in params]))
    nfm += '> %d' %  (sum([ int(param) for param in params]))
    mixed = []
    mixed.append(nfm)
    mixed.append('?' + 'getspar' + ':0=%d> 0' %random.randrange(0,5))

    return response_list('getspar', cmd_num, app_num, expected_ack, mixed, response_type)


def response_list(cmd, cmd_num, app_num, expected_ack, mixed, response_type):
    """This function make an answer for a request.

    Parameters:

    - `cmd`: the command (getpos, getspar, etc.)
    - `cmd_num`: the command identification number
    - `app_num`: the application number (minor servo address)
    - `expected_ack`: the expected_ack response
    - `mixed`: a mixed answer
    - `response_type`: the type of response: expected_ack, mixed

    Return a list of answers.
    """
    answers = []
    if response_type == "expected_ack":
        answers.append(expected_ack + closers[0])
    elif response_type == "mixed":
        for message in mixed:
            answers.append(message + closers[0])
        if(random.randrange(0,10) == 5):
            answers.append(('!NAK_' + cmd + ':%d=%d> Wrong State = Application ...' %(cmd_num, app_num)) + closers[0])
        else:
            answers.append(('!NAK_' + cmd + ':0=%d> Wrong State = Application ...' %app_num) + closers[0])
        answers.append('invalid')
        answers.append('?invalid')
        answers.append(expected_ack + closers[0])
        answers.append('?invalid' + closers[0])

    return answers