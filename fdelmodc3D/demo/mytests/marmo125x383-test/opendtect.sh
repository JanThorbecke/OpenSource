#!/bin/bash

ODTPATH=/home/victork/OpendTect/6.4.0

here=$(pwd)

$ODTPATH/mk_datadir $here
$ODTPATH/start_dtect 
