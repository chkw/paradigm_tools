#!/usr/bin/awk -f
BEGIN{FS="\t";}
{if (NF < 3) {print;};}
END{}
