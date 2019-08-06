#!/bin/bash

for e in GW150914 GW151226 GW170104 GW170608 GW170729 GW170809 GW170814 \
         GW170817 GW170818 GW170823 S190412m S190421ar S190425z S190426c \
         S190503bf S190510g S190512at S190513bm S190517h S190519bj S190521g \
         S190521r; do 
  if ! [ -e ../ligobgresults-$e/fardet-t02.list ]; then
    echo ../ligobgresults-$e/fardet-t02.list does not exist, skipping
    continue
  fi

  if ! [ -e  ../ligosidebandresults-$e/*-fardet-t02.hadded.root ]; then
    echo ../ligosidebandresults-$e/*-fardet-t02.hadded.root does not exist, skipping
    continue
  fi

  ./tricia.sh ../ligobgresults-$e/fardet-t02.list \
    ../ligosidebandresults-$e/*-fardet-t02.hadded.root
done

# note
exit 0


for e in GW150914 GW151226 GW170104 GW170608 GW170729 GW170809 \
         GW170814 GW170818 S190421ar S190426c S190513bm S190521g S190521r; do
  ./tricia.sh ../ligobgresults-$e/neardet-ddactivity.list \
    ../ligosidebandresults-$e/*-neardet-ddactivity1.hadded.root
done

for f in \
ligosidebandresults-GW151226/2015-12-26T02:38:53.65Z-fardet-ddenergy.hadded.root \
ligosidebandresults-GW170104/2017-01-04T09:11:58.6Z-fardet-ddenergy.hadded.root \
ligosidebandresults-GW170608/2017-06-08T01:01:16.49Z-fardet-ddenergy.hadded.root \
ligosidebandresults-GW170729/2017-07-29T17:56:29.3Z-fardet-ddenergy.hadded.root \
ligosidebandresults-GW170809/2017-08-09T07:28:21.8Z-fardet-ddenergy.hadded.root \
ligosidebandresults-GW170814/2017-08-14T09:30:43.53Z-fardet-ddenergy.hadded.root \
ligosidebandresults-GW170817/2017-08-17T11:41:04.4Z-fardet-ddenergy.hadded.root \
ligosidebandresults-GW170818/2017-08-18T01:25:09.1Z-fardet-ddenergy.hadded.root \
ligosidebandresults-GW170823/2017-08-23T12:13:58.5Z-fardet-ddenergy.hadded.root \
ligosidebandresults-S190412m/2019-04-12T04:30:44.165622Z-fardet-ddenergy.hadded.root \
ligosidebandresults-S190421ar/2019-04-21T20:38:56.250977Z-fardet-ddenergy.hadded.root \
ligosidebandresults-S190425z/2019-04-25T07:18:05.017147Z-fardet-ddenergy.hadded.root \
ligosidebandresults-S190426c/2019-04-26T14:21:55.336540Z-fardet-ddenergy.hadded.root \
ligosidebandresults-S190503bf/2019-05-03T17:54:04.294490Z-fardet-ddenergy.hadded.root \
ligosidebandresults-S190510g/2019-05-10T01:59:39.291636Z-fardet-ddenergy.hadded.root \
ligosidebandresults-S190512at/2019-05-12T17:07:14.422363Z-fardet-ddenergy.hadded.root \
ligosidebandresults-S190513bm/2019-05-13T19:54:28.747089Z-fardet-ddenergy.hadded.root \
ligosidebandresults-S190517h/2019-05-17T04:51:01.830582Z-fardet-ddenergy.hadded.root \
ligosidebandresults-S190519bj/2019-05-19T14:35:44.397949Z-fardet-ddenergy.hadded.root \
ligosidebandresults-S190521g/2019-05-21T02:02:29.447266Z-fardet-ddenergy.hadded.root \
ligosidebandresults-S190521r/2019-05-21T06:43:59.463379Z-fardet-ddenergy.hadded.root; do
  ./tricia.sh ../ligobgresults/fardet-ddenergy.list ../$f
done
