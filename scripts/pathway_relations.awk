#!/usr/bin/awk -f
BEGIN{FS="\t";}
{if (NF > 2) {print;};}
END{}
