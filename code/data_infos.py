



cell_types = ["MSCH", "THP1", "MSCS", "Fibro", "MSCSmapH"]

table_columns = {"THP1": {
                            "2D": ['LFQI 2D '+str(i) for i in range(1,4)],
                            "3D": ['LFQI 3D '+str(i) for i in range(1,4)],
                            "MS": ['LFQI MS '+str(i) for i in range(1,4)],
                            "HS": ['LFQI HS '+str(i) for i in range(1,4)],
                            "MS.6h": ['LFQI MS.6h '+str(i) for i in range(1,4)],
                            
                            "Lysat MS": ["LFQI Lysat MS "+str(i) for i in range(1,4)],
                            "Lysat MS.6h": ["LFQI Lysat MS.6h "+str(i) for i in range(1,4)],
                            "Lysat HS": ["LFQI Lysat HS "+str(i) for i in range(1,4)],
                            "Lysat 3D": ["LFQI Lysat 3D "+str(i) for i in range(1,4)],
                            "Lysat 2D": ["LFQI Lysat 2D "+str(i) for i in range(1,4)],
                         
                            },
            "MSCH": {
                            "2D": ['LFQI 2D '+str(i) for i in range(1,4)],
                            "3D": ['LFQI 3D '+str(i) for i in range(1,4)],
                            "HS": ['LFQI HS '+str(i) for i in range(1,4)],
                    },
            "MSCS": {
                            "2D": ['LFQI 2D '+str(i) for i in range(1,4)],
                            "3D": ['LFQI 3D '+str(i) for i in range(1,4)],
                            "MS": ['LFQI MS '+str(i) for i in range(1,4)],                            
                            "HS": ['LFQI HS '+str(i) for i in range(1,4)],
                    },
            "Fibro": {
                            "2D": ['LFQI 2D '+str(i) for i in range(1,4)],
                            "3D": ['LFQI 3D '+str(i) for i in range(1,4)],
                            "HS": ['LFQI HS '+str(i) for i in range(1,4)],
                    },
                    

            }
table_columns["MSCSmapH"] = table_columns["MSCS"]


protocols = {"THP1": ["2D","3D","MS","HS", "MS.6h"],
            "MSCH": ["2D","3D","HS"],
            "MSCS": ["2D","3D","MS","HS"],
            "Fibro": ["2D","3D","HS"],
            }
protocols["MSCSmapH"] = protocols["MSCS"]


protocols_lysat = {"THP1": ["Lysat 2D", "Lysat 3D","Lysat MS","Lysat HS", "Lysat MS.6h"],
            "MSCH": [],
            "MSCS": [],
            "Fibro": [],
            }
protocols_lysat["MSCSmapH"] = protocols_lysat["MSCS"]


experiments = {cell_type: sum([table_columns[cell_type][protocol] \
                               for protocol in protocols[cell_type]],[]) \
                                for cell_type in cell_types}

experiments_lysat = {cell_type: sum([table_columns[cell_type][protocol] \
                               for protocol in protocols_lysat[cell_type]],[]) \
                                for cell_type in cell_types}


all_protocols = sum([[k+"_"+protocol  for protocol in v] for k, v in protocols.items()],[])




dic_species = {"MSCH":"homo_sapiens",
               "THP1":"homo_sapiens",
               "Fibro":"homo_sapiens",
               "MSCS":"mus_musculus",
               "MSCSmapH": "mus_musculus"}

dic_protocol_corresp = {"THP1":{},
                        "MSCH":{},
                        "Fibro":{},
                        "MSCS":{}}

        
for measure in range(1,4):
       
    dic_protocol_corresp["THP1"]["LFQ intensity EV-300RPM-3h_"+str(measure)] = "LFQI MS "+str(measure)
    dic_protocol_corresp["THP1"]["LFQ intensity EV-300RPM-6h_"+str(measure)] = "LFQI MS.6h "+str(measure)
    dic_protocol_corresp["THP1"]["LFQ intensity EV-500RPM-3h_"+str(measure)] = "LFQI HS "+str(measure)
    dic_protocol_corresp["THP1"]["LFQ intensity EV-car3D-72h_"+str(measure)] = "LFQI 3D "+str(measure)
        
dic_protocol_corresp["THP1"]["LFQ intensity EV-car2D-72h_3"] = "LFQI 2D 3"
dic_protocol_corresp["THP1"]["LFQ intensity EV-car2D-72hDR_1"] = "LFQI 2D 1"
dic_protocol_corresp["THP1"]["LFQ intensity EV-car2D-72hDR_2"] = "LFQI 2D 2"               
                

for measure in range(1,4):
    dic_protocol_corresp["THP1"]["LFQ intensity L-300RPM-3h_"+str(measure)] = "LFQI Lysat MS "+str(measure)
    dic_protocol_corresp["THP1"]["LFQ intensity L-300RPM-6h_"+str(measure)] = "LFQI Lysat MS.6h "+str(measure)
    dic_protocol_corresp["THP1"]["LFQ intensity L-500RPM-3h_"+str(measure)] = "LFQI Lysat HS "+str(measure)
    dic_protocol_corresp["THP1"]["LFQ intensity L-car3D-72h_"+str(measure)] = "LFQI Lysat 3D "+str(measure)
    dic_protocol_corresp["THP1"]["LFQ intensity L-car2D-72h_"+str(measure)] = "LFQI Lysat 2D "+str(measure)


for measure in range(1,4):
    dic_protocol_corresp["MSCS"]["LFQ intensity 2D"+str(measure)] = "LFQI 2D "+str(measure)
    dic_protocol_corresp["MSCS"]["LFQ intensity 3D"+str(measure)] = "LFQI 3D "+str(measure)
    dic_protocol_corresp["MSCS"]["LFQ intensity HS"+str(measure)] = "LFQI HS "+str(measure)
    dic_protocol_corresp["MSCS"]["LFQ intensity MS"+str(measure)] = "LFQI MS "+str(measure)


for measure in range(1,4):
    dic_protocol_corresp["MSCH"]["LFQ intensity 2D_"+str(measure)] = "LFQI 2D "+str(measure)
    dic_protocol_corresp["MSCH"]["LFQ intensity 3D_4"+str(measure)] = "LFQI 3D "+str(measure)
    dic_protocol_corresp["MSCH"]["LFQ intensity HighSpeed_"+str(37 + measure)] = "LFQI HS "+str(measure)



for measure in range(1,4):
    dic_protocol_corresp["Fibro"]["LFQ intensity EV_3D"+str(measure)] = "LFQI 3D "+str(measure)
    dic_protocol_corresp["Fibro"]["LFQ intensity EV_HS"+str(measure)] = "LFQI HS "+str(measure)

dic_protocol_corresp["Fibro"]["LFQ intensity EV_2D1"] = "LFQI 2D 1"
dic_protocol_corresp["Fibro"]["LFQ intensity EV_2D2DR"] = "LFQI 2D 2"
dic_protocol_corresp["Fibro"]["LFQ intensity EV_2D3"] = "LFQI 2D 3"




# list_comparisons = [("2D", "3D"), ("2D", "HS"), ("3D", "HS"), ("2D", "MS"), ("3D", "MS"), ("MS", "HS"), ("2D", "MS.6h"), ("3D", "MS.6h"), ("MS", "MS.6h"), ("MS.6h", "HS")]



# interesting_comparisons = {cell_type: [comparison for comparison in list_comparisons if comparison[0] in protocols[cell_type] and comparison[1] in protocols[cell_type]] for cell_type in cell_types}



dic_protocol_corresp["MSCSmapH"] = dic_protocol_corresp["MSCS"]
experiments["MSCSmapH"] = experiments["MSCS"]

        