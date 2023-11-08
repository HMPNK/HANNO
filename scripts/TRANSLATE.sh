#!/bin/bash
cat $1 | perl $SCRIPTS/translate.perl | grep -v ^$ 
