# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 16:49:20 2021

@author: Matteo
"""

import os

def ipFinder(gateway, start_search, stop_search):
    #Gateway is of type 192.168.1.1
    str_gate = gateway.split('.')
    base_id = str_gate[:-1]
    new_id = str()
    for ele in base_id:
        new_id += ele+'.'
    
    for num in range(start_search, stop_search):
        ip_to_check = new_id + str(num)
        print(ip_to_check)
        result = os.system('ping -n 1 '+ ip_to_check)
        if result == 0:
            print(ip_to_check + ' ' + 'Found!')
            
            

if __name__ == '__main__':
    ipFinder('192.168.1.1',100,200)
    