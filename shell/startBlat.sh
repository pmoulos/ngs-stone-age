#!/bin/bash

# Name this abomination
SCRIPT_NAME=$(basename $0)

# We must be under sudo
if [ `id -u` -ne 0 ]; 
then
	echo "You must execute this script as root!"
	exit 1
fi

# Defaults
START=$1
GBDB_PATH=/media/HD2/gbdb
KENTHOME=/media/HD2/kenthome/bin/x86_64
HOST=localhost
PORT=60000
DBUSER=gbuser

if [ $START != "first" -a $START != "start"  -a $START != "stop" ]
then
	echo ""
	echo "The argument must be one of first, start or stop!"
	echo ""
	exit 1
fi

do_sql()
{
	echo ""
	TMPORT=$3
	MYSQLCMD="mysql -u"$1" -p -e \"use hgcentral; "
	for ORG in `dir -d *`
	do
		if [ $ORG != "hgFixed" -a $ORG != "genbank" ]
		then
			MYSQLCMD=$MYSQLCMD"UPDATE blatServers SET host='"$2"', port="$TMPORT" WHERE db LIKE '%"$ORG"%'; "
			TMPORT=$(($TMPORT+1))
		fi
	done
	MYSQLCMD=$MYSQLCMD"\""
	echo "Updating hgcentral blat servers... You might be asked for a password"
	echo "$MYSQLCMD"
	#$MYSQLCMD
}

do_blat()
{
	echo ""
	TMPORT=$3
	for ORG in `dir $GBDB_PATH`
	do
	if [ $ORG != "hgFixed" -a $ORG != "genbank" ]
		then
			cd $GBDB_PATH/$ORG
			echo "Starting blat server for $ORG... Please wait..."			
			#echo "$1/gfServer start $2 $TMPORT -stepSize=5 -log=/tmp/blatlog.$ORG.log $ORG.2bit"
			$1/gfServer start $HOST $TMPORT -stepSize=5 -log=/tmp/blatlog.$ORG.log $ORG.2bit &
			TMPORT=$(($TMPORT+1))
		fi
	done
}

kill_blat()
{
	for PID in `pidof gfServer`
	do
		kill -9 $PID
	done
}

# The actual work
if [ $START = "first" ]
then
	do_sql $DBUSER $HOST $PORT
	do_blat $KENTHOME $HOST $PORT
elif [ $START = "start" ]
then
	do_blat $KENTHOME $HOST $PORT
elif [ $START = "stop" ]
	kill_blat
fi
echo ""
