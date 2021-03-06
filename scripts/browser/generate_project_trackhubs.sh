# !/bin/bash

# this is set to only grab the first 1000 projects and aggregations
# but there shouldn't be more than 1000 active projects and the browser is impractical with 1000 aggregations

# in the future, only remake things that have changed jsons

configfile=/home/solexa/stampipes-hpc/config/ucsc_browser/trackhub_projects.config
priority=400 # legacy variable

cd /net/seq/data/aggregations/project_trackhubs/

# pull all projects
python $STAMPIPES/scripts/browser/parse_all_projects.py -o project_list.txt

# for each project
while read line; do
	id=$(echo "$line" | cut -f 1)
	name=$(echo "$line" | cut -f 2)
	name=$(echo $name | tr ' ' '_' | tr '/' '-')
	python $STAMPIPES/scripts/lims/get_processing.py -p $id -o project_${id}_${name}.json
	python $STAMPIPES/scripts/browser/make_trackhubs_for_projects.py -c $configfile -p $priority -j project_${id}_${name}.json -n ${name}
done < project_list.txt
