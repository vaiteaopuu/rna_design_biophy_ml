rfam = ["RF00005_1" , "RF01852" , "RF00059" , "RF00100" , "RF00162" , "RF00167"
, "RF00168" , "RF00234" , "RF00380" , "RF01051" , "RF01510" , "RF01734" ,
"RF01750" , "RF01763" , "RF01786" , "RF01831_1" , "RF01854" , "RF00028_1" ,
"RF00169" , "RF01767" , "RF02682" , "RF00163" , "RF01054" , "RF02695" ,
"RF00233" , "RF01807" , "RF01415" , "RF00379" , "RF02683" , "RF00921" ,
"RF02679" , "RF01982" , "RF00606_1" , "RF01300" , "RF02681" , "RF00442_1" ,
"RF00080" , "RF01826" , "RF03017" , "RF00011_1" , "RF02540" , "RF00029" ,
"RF00010" , "RF00164" , "RF00458" , "RF00050" , "RF00504" , "RF00044" ,
"RF00061_2" , "RF01689" , "RF02447" , "RF01725" , "RF02266" , "RF02927" ,
"RF02680" , "RF02001_2" , "RF02553"]

pdb = ["1ehz_A" , "3rg5_A" , "3d2g_A" , "5lys_A" , "3gx5_A" , "4tzx_X" ,
        "3dil_A" , "2h0s_B" , "3pdr_A" , "4yaz_A" , "3slq_A" , "4enc_A" , "4xwf_A" ,
        "5nwq_A" , "3q3z_A" , "4lvv_A" , "4wfl_A" , "1gid_A" , "1z43_A" , "3e5c_A" ,
        "3nkb_B" , "3zp8_A" , "4jf2_A" , "4k27_U" , "4p5j_A" , "4p95_A" , "4pqv_A" ,
        "4qln_A" , "4rum_A" , "5dun_A" , "5k7d_A" , "5kpy_A" , "5m0h_A" , "5ob3_A" ,
        "5t5a_A" , "5u3g_B" , "6cb3_A" , "6fz0_A" , "2oiu" , "1nbs_A" , "1ffz_A" ,
        "1kxk" , "1u9s" , "1xjr" , "2il9" , "3f2q" , "3ox0_A" , "3r4f" , "3t4b_A" ,
        "4frg" , "4jrc_A" , "4l81" , "4plx_A" , "4r4v" , "4rzd_A" , "4y1o_A" , "6cu1_A"]

map_rfam_pdb = {fam: struct for fam, struct in zip(rfam, pdb)}

rfam_test = ["RF00001" ,"RF00008" ,"RF00017" ,"RF00023" ,"RF00027" ,"RF00102" ,"RF00166" ,"RF00174" ,"RF00207" ,"RF00209" ,"RF00210" ,"RF00374" ,"RF00480" ,"RF00500" ,"RF00634" ,"RF01073" ,"RF01381" ,"RF01704" ,"RF01727" ,"RF01739" ,"RF01857" ,"RF01998" ,"RF02012"]
pdb_test = ["1c2x" ,"5di2" ,"1l9a" ,"1p6v" ,"5zal" ,"6ol3" ,"2mf0" ,"4gma" ,"2ke6" ,"4c4q" ,"2nbx" ,"1s9s" ,"1z2j" ,"2krl" ,"6ues" ,"2lc8" ,"2n1q" ,"6qn3" ,"6hag" ,"5ddp" ,"3ndb" ,"3bwp" ,"4r0d"]
map_rfam_pdb_test = {fam: struct for fam, struct in zip(rfam_test, pdb_test)}

def read_wt_seq(infile):
    "read from coconet dataset"
    seq = ""
    for l in open(infile):
        l.replace(" ", "")
        if l.startswith("SEQUENCE"):
            seq += l.strip().split(":")[1].strip().upper()
    return seq
