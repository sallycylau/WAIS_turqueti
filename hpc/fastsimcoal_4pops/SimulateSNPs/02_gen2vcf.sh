cat 1_WS_AS_RS_EA_psccc_fulcol1.out/1_WS_AS_RS_EA_psccc_fulcol1.out_1_1.gen | awk -f gen2vcf.awk > 1_WS_AS_RS_EA_psccc_fulcol1.out.vcf

cat 1_WS_AS_RS_EA_psc_nocol.out/1_WS_AS_RS_EA_psc_nocol.out_1_1.gen | awk -f gen2vcf.awk > 1_WS_AS_RS_EA_psc_nocol.out.vcf

cat 1_WS_AS_RS_EA_psc_parcol.out/1_WS_AS_RS_EA_psc_parcol.out_1_1.gen | awk -f gen2vcf.awk > 1_WS_AS_RS_EA_psc_parcol.out.vcf
