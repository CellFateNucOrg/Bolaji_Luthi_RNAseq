tea -d ../.././publicData/anatomy_dict_95_33_WS285.csv -q 0.05 -s COH1cs_upGenes_WBID.txt COH1cs_up_tissue tissue
tea -d ../.././publicData/phenotype_dict_95_50_WS285.csv -q 0.05 -s COH1cs_upGenes_WBID.txt COH1cs_up_phe phenotype
tea -d ../.././publicData/go_dict_95_33_WS285.csv -q 0.05 -s COH1cs_upGenes_WBID.txt COH1cs_up_go go
tea -d ../.././publicData/anatomy_dict_95_33_WS285.csv -q 0.05 -s COH1cs_downGenes_WBID.txt COH1cs_down_tissue tissue
tea -d ../.././publicData/phenotype_dict_95_50_WS285.csv -q 0.05 -s COH1cs_downGenes_WBID.txt COH1cs_down_phe phenotype
tea -d ../.././publicData/go_dict_95_33_WS285.csv -q 0.05 -s COH1cs_downGenes_WBID.txt COH1cs_down_go go
cd ../../
