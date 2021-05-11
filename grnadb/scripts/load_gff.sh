#!/bin/sh

#/xraid/bioware/linux/bioperl_scripts/bp_bulk_load_gff.pl --create --user=gid --password gidgid123 --fasta gff --local --database dbi:mysql:$1:gmoddb.mbl.edu gff/*
/xraid/habitat/mcipriano/perllib/bin/bp_bulk_load_gff.pl --create --user=gid --password gidgid123 --fasta gff --local --database dbi:mysql:$1:gmoddb.mbl.edu gff/*
