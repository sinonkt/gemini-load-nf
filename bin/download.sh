for db in GRCh37.75 GRCm38.75; do 
  java -jar $SNPEFF_JAR download $db; 
done