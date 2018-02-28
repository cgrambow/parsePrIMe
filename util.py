#!/usr/bin/env python
# -*- coding:utf-8 -*-

import cPickle as pickle


def custom_warning(msg, *args, **kwargs):
    return 'Warning: ' + str(msg) + '\n'


def pickle_dump(path, obj):
    with open(path, 'wb') as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)


def pickle_load(path):
    with open(path, 'rb') as f:
        return pickle.load(f)
