for file in pca_table_split/*; do echo "bsubl -oo filter.merge."${file##*pca_table_split/}".log 'bash filter_merge_germline_by_sample.sh "$file"'"; done
