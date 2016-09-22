#!/usr/bin/python

import os

os.system("nohup ./m3birfinder config­1.txt > ./Log/run0222.log 2>&1")
os.system("./create_columns.sh")
os.system("mysql pcilnoothbec5 mmbirfinder")

print("Execute these commands in mysql:")
print("----------------------------------------------------------------------")
print("delete from mmbirfinder.unaligned;")
print("LOAD DATA LOCAL INFILE ‘columns.txt’ into table mmbirfinder.unaligned;")
print("exit")
print("----------------------------------------------------------------------")

os.system("sudo -u mysql /bin/myisamchk -o /var/lib/mysql/mmbirfinder/unaligned.MYI --key_buffer_size=2G")

os.system("cd trial2")
os.system("perl ../db_fixed.pl")
os.system("cd ..")
os.system("ln -s trial2/mysql_results.txt .")
os.system("./m3birfinder config-2.txt > ./Log/run0222_1.log 2>&1")
