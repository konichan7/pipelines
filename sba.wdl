#CAPER docker quay.io/encode-dcc/atac-seq-pipeline:v1.4.2

workflow cross_corr {

    Array[File] rna_bws
    Array[File] atac_bws
    Array[Pair[File, File]] pairs = cross(rna_bws, atac_bws)
    
    scatter(pair in pairs){
        File rna_bw = pair.left
        File atac_bw = pair.right
        call sba1 { input: in1=rna_bw, in2=atac_bw }
        call sba2 { input: in1=sba1.out1, in2=sba1.out2 }
    }
}

task sba1{
    
    File in1
    File in2
    String chr_atac
    String chr_rna

    command <<<
        python <<CODE
        
        #!/usr/bin/env python
        # coding: utf-8
        
        import numpy
        import pandas
        import pyBigWig
        import matplotlib.pyplot as plt

        rna_bw = pyBigWig.open("${in1}") #bigwig.file
        atac_bw = pyBigWig.open("${in2}") #bigwig.file

        rna_intervals = rna_bw.intervals('${chr_rna}')
        atac_intervals = atac_bw.intervals('${chr_atac}')
        out1 = os.path.basename(rna_bw).strip('.bigWig') + '.rna.qc.csv'
        out2 = os.path.basename(atac_bw).strip('.bigWig') + '.atac.qc.csv'

        mean_s=rna_bw.stats('${chr_rna}')[0]
        max_s=rna_bw.stats('${chr_rna}',type="max",nBins=3)

        peak_all=[]

        for interval in rna_intervals:
            if interval[2] in max_s:
                peaks=[]
                for i in range(interval[0]+1,interval[1]+1):
                    peaks.append((i,interval[2]))
                peak_all.append(peaks)

        rna_peaks=[]

        for regions in peak_all:
            front = [(i,rna_bw.stats('${chr_rna}',i-1,i)[0]) for i in range(regions[0][0]-500, regions[0][0])] #chr19
            end = [(i,rna_bw.stats('${chr_rna}',i-1,i)[0]) for j in range(regions[-1][0]+1, regions[-1][0]+501)] #chr19
            rna_peaks.append(front+regions+end)

        means = atac_bw.stats('${chr_atac}')[0] #chr19
        maxs = atac_bw.stats('${chr_atac}',type="max",nBins=22) #chr19

        all_peak=[]

        for interval2 in atac_intervals:
            if interval2[2] in maxs:
                peaks=[]
                for i in range(interval2[0]+1,interval2[1]+1):
                    peaks.append((i,interval2[2]))
                all_peak.append(peaks)

        atac_peaks=[]

        for regions in all_peak:
            front=[(i,atac_bw.stats('${chr_atac}',i-1,i)[0]) for i in range(regions[0][0]-500,regions[0][0])] #chr19
            end=[(i,atac_bw.stats('${chr_atac}',i-1,i)[0]) for i in range(regions[-1][0]+1,regions[-1][0]+501)] #chr19
            atac_peaks.append(front+regions+end)

        qc_atac = pandas.DataFrame(atac_peaks)
        qc_rna = pandas.DataFrame(rna_peaks)

        qc_atac.to_csv(out1,index=False) #csv File name
        qc_rna.to_csv(out2,index=False) #csv File name
        
        CODE
    >>>
    output{
        File out1 = glob("*.rna.qc.csv")[0]
        File out2 = glob("*.atac.qc.csv")[0]
    }
    runtime{
        cpu:2
        memory:"8000MB"
        disks : "/mnt/disks/sydir/ 40 SSD"
        docker: "quay.io/encode-dcc/atac-seq-pipeline:v1.4.2"
    }
}

task sba2{
    
    File in1
    File in2

    command <<<
        python<<CODE
    
        #!/usr/bin/env python
        # coding: utf-8

        import inspect as ins
        import numpy as np
        import pandas as pd
        import pyBigWig
        import pysam
        import matplotlib.pyplot as plt
        import scipy.stats as st


        rna_seq = pd.read_csv("${in1}")
        atac_seq = pd.read_csv("${in2}")

        rna=[]
        whole_qc_rna=[]
        atac_qc=[]
        whole_qc_atac=[]

        for i in range(rna_seq.index.stop):
            rna.append(list(rna_seq.iloc[i,:].dropna()))

        def str_to_ls(input_str):
            str_list = input_str[1:-1].split(',')
            new_list = []
            new_list.append(int(str_list[0]))
            if str_list[1] != ' None':
                new_list.append(float(str_list[1]))
            else:
                new_list.append(0)
            return(new_list)

        for j in range(len(rna)):
            ls=[]
            for i in range(len(rna[j])):
                ls.append(str_to_ls(rna[j][i]))
            whole_qc_rna.append(ls)

        for i in range(atac_seq.shape[0]):
            atac_qc.append(list(atac_seq.iloc[i,:].dropna()))

        for j in range(len(atac_qc)):
            ls = []
            for i in range(len(atac_qc[j])):
                ls.append(str_to_ls(atac_qc[j][i]))
            whole_qc_atac.append(ls)

        def where_skewed(sig_list):
            sig_list = pd.DataFrame(sig_list)
            x = sig_list.values[:,0]
            y = sig_list.values[:,1]
            if st.skew(y) > 0:
                return(1)                   # right skewed
            elif st.skew(y) == 0:
                return(0)                   # symmetric
            else:
                return(-1)                  # left skewed

        def len_saming(rna_list, atac_list):
            new_list=[]
            if len(rna_list) > len(atac_list):
                diff = len(rna_list) - len(atac_list)
                atac_list = pd.DataFrame(atac_list)
                x = atac_list.values[:,0]
                y = atac_list.values[:,1]
                skewness = st.skew(y)
                if where_skewed(atac_list) == 1:                          #right skewed -> skew : 1
                    add_pre = [0]*int(np.around(diff*skewness/(1+skewness)))
                    add_app = [0]*(diff - int(np.around(diff*skewness/(1+skewness))))

                elif where_skewed(atac_list) == -1:                                     #left skewed -> 1 : -skew
                    add_pre = [0]*int(np.around(diff/(1-skewness)))
                    add_app = [0]*(diff - int(np.around(diff/(1-skewness))))

                else:
                    add_pre = [0]*int(np.around(diff/2))
                    add_app = [0]*(diff - int(np.around(diff/2)))

                y = np.insert(y,0,add_pre)
                y = np.append(y,add_app)
                rna_list = pd.DataFrame(rna_list)
                y2 = rna_list.values[:,1]
                return(y2,y)
            else:
                diff = len(atac_list) - len(rna_list)
                rna_list = pd.DataFrame(rna_list)
                x = rna_list.values[:,0]
                y = rna_list.values[:,1]
                skewness = st.skew(y)
                if where_skewed(rna_list) == 1:                              #right skewed -> skew : 1
                    add_pre = [0]*int(np.around(diff*skewness/(1+skewness)))
                    add_app = [0]*(diff - int(np.around(diff*skewness/(1+skewness))))

                elif where_skewed(rna_list) == -1:                           #left skewed -> 1 : -skew
                    add_pre = [0]*int(np.around(diff/(1-skewness)))
                    add_app = [0]*(diff - int(np.around(diff/(1-skewness))))

                else:
                    add_pre = [0]*int(np.around(diff/2))
                    add_app = [0]*(diff - int(np.around(diff/2)))

                y = np.insert(y,0,add_pre)
                y = np.append(y,add_app)
                atac_list = pd.DataFrame(atac_list)
                y2 = atac_list.values[:,1]
                return(y,y2)

        def cal_corr(rna,atac):
            corr_mat = []
            index = []
            columns = []
            for i in range(len(atac)):
                index.append("atac[%d]" %(i))

            for i in range(len(rna)):
                columns.append("rna[%d]" %(i))

            for atac_region in atac:
                corr_row = []
                for rna_region in rna:
                    if len(rna_region) == len(atac_region):
                        rna_region = pd.DataFrame(rna_region)
                        atac_region = pd.DataFrame(atac_region)
                        yr = rna_region.values[:,1]
                        ya = atac_region.values[:,1]
                        corr, p_val = st.pearsonr(yr, ya)
                        #if np.isnan(corr):
                        #    corr = 0
                        #    p_val = 1
                        #corr_row.append(corr)
                    else:
                        corr, p_val = st.pearsonr(len_saming(rna_region,atac_region)[0],len_saming(rna_region,atac_region)[1])
                        #if np.isnan(corr):
                        #    corr = 0
                        #    p_val = 1
                        #corr_row.append(corr)
                    corr_row.append(corr)
                corr_mat.append(corr_row)
            corr_mat = pd.DataFrame(corr_mat, index = index, columns = columns)
            return(corr_mat)

        result = cal_corr(whole_qc_rna,whole_qc_atac)

        result.to_csv("result_test.csv")

        CODE
    >>>
    output{
        File out = "result_test.csv"
    }
    
    runtime{
        cpu : 2
        memory : "8000MB"
        disks : "/mnt/disks/sydir 40 SSD"
        docker: "quay.io/encode-dcc/atac-seq-pipeline:v1.4.2"
    }
}

