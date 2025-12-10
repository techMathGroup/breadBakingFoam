#!/usr/bin/python

# FILE DESCRIPTION======================================================

# Python script to create docs folder from Github wiki for Github pages

# IMPORTS===============================================================
import os
import shutil
import re

# SETUP=================================================================
wikiGithub = "https://github.com/hlavatytomas/breadBakingFoam.wiki.git"

# CODE==================================================================
# -- get the name of the wiki page
nameOfWiki = wikiGithub.split('/')[-1].replace('.git', '')
# -- clone the repository with the wiki
if os.path.exists(nameOfWiki):
    shutil.rmtree(nameOfWiki)
os.system('git clone ' + wikiGithub)

# -- go through all the files in wiki, copy them in docs and correct 
# -- the paths
files = [file.name for file in os.scandir(nameOfWiki) if file.is_file() ]
for file in files:
    shutil.copyfile(nameOfWiki + '/' + file, file)

    with open(file, "r", encoding="utf8") as f:
        content = f.read()

    content = re.sub(r'(src=")[^"]*/docs/([^"]+)(")', r'\1\2\3', content)

    with open(file, "w", encoding="utf8") as f:
        f.write(content)

# -- remove wiki clone
shutil.rmtree(nameOfWiki)