# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 14:11:26 2022

@author: ellio
"""

import sqlite3
import numpy as np
import matplotlib.pyplot as plt

t=0
u=1
v=2
p=3

conn = sqlite3.connect('données_mesure.db')

cursor = conn.cursor()

cursor.execute("""
                   INSERT INTO valeurs(condition ,horizontal ,vertical ,pressure) VALUES (?,?,?,?)""",(t,u,v,p))

c=cursor.execute("""
                 SELECT condition,horizontal,vertical,pressure from valeurs""")

for row in c:
    for i in range(0,3):
        print(row[i])

conn.commit()