import matplotlib.pyplot as plt
import pandas as pd
import pysam
import os
import ast
import numpy as np
import json
import collections
from tqdm import tqdm
from Bio.Seq import Seq
from Bio import SeqIO
import more_itertools as mit
import argparse as ap


def find_sequence(sequence, subsequence):
    positions = []
    start = 0
    while True:
        start = sequence.find(subsequence, start)
        if start == -1:
            break
        positions.append(start)
        start += 1
    return positions




def orientationFinder(df):
    df2 = df.copy()
    df2['Orientation']='NONE'
    for row in df2.index:
        #print(row)
        if df2.at[row,'TE_Hits'] == 'NONE':
            #print("NO HITS")
            continue
        else:
            elementAnnotation = str(df2.at[row,'Element_Annotation'])
            designation = str(df2.at[row,'TE_Designation'])
            sense=0
            antisense=0
            hitList = ast.literal_eval(str(df2.at[row,'TE_Hits']))
            senseSize =0
            antisenseSize=0
            
            count = 0
            umbrellaFlag=0
            for hit in hitList:
                umbrella= str(hit.split()[10])
                
                if str(hit.split()[9]) == elementAnnotation:
                    count+=1
                    
                    if str(hit.split()[8])=='+':
                        sense+=1
                        senseSize+=int(hit.split()[6])-int(hit.split()[5])
                    else:
                        antisenseSize+=int(hit.split()[6])-int(hit.split()[5])
                        antisense+=1
                        
                elif count ==0 and (umbrella == 'Simple_repeat' or umbrella == 'Low_complexity'):
                    count+=1
                    antisense+=1
                    umbrellaFlag=1
                
                elif count == (len(hitList)-1) and (umbrella == 'Simple_repeat' or umbrella == 'Low_complexity'):
                    count+=1
                    sense+=1
                    umbrellaFlag=1
                    
                elif designation == umbrella:
                    count+=1
                    if str(hit.split()[8])=='+':
                        sense+=.5
                        senseSize+=int(hit.split()[6])-int(hit.split()[5])
                    else:
                        antisenseSize+=int(hit.split()[6])-int(hit.split()[5])
                        antisense+=.5
                
                else:
                    count+=1
                    continue
                    
            if umbrellaFlag>0:
                pass
            else:
                if senseSize>=antisenseSize:
                    sense+=1
                else:
                    antisense+=1
            
            if sense>0 and antisense ==0:
                df2.at[row,'Orientation']='+'
            elif sense==0 and antisense>0:
                df2.at[row,'Orientation']='-'
            elif sense>antisense:
                df2.at[row,'Orientation']='+'
            elif antisense>sense:
                df2.at[row,'Orientation']='-'
            else:
                continue
    return(df2)

    


def tailCounter(df):
    
    df2= df.copy()
    df2['FILTER_RESULTS']='Good_Row'
    df2['Tail_Begins']='No_Tail_Detected'
    df2['Tail_Type']='No_Tail_Type'
    df2['Tail_Length']=0
    df2['Tail_Seed_Hits']=0
    
    for row in df2.index:
        

        sequence = str(df2.at[row,'Sequence']).upper()
        
        tTail = 'TTTTT'
        tTailFlag=0
        trunningTotal = 0
        tTailList=[]
        
        if tTail in sequence:
            taillength=0
            findings = [int(x) for x in find_sequence(sequence, tTail)]
            if len(findings)>0:
                tTailFlag=1
                tmin = min(findings)
                
                iterable = findings
                df2.at[row,'Tail_Seed_Hits']=len(findings)
                groupings = [list(group) for group in mit.consecutive_groups(iterable)]
                flag=0
                for group in groupings:
                    if flag==0:

                        taillength+= len(tTail) + (len(group)-1)
                        groupEnd = max(group)+len(tTail)
                        flag=1

                    else:

                        if min(group)-groupEnd<=4:
                            taillength+= (len(tTail) + (len(group)-1) + abs(min(group)-groupEnd))
                            groupEnd = max(group)+len(tTail)
                        else:
                            continue

                tTailLength = taillength
                
            else:
                pass
                
        
        aTail = 'AAAAA'
        aTailFlag=0
        arunningTotal = 0
        aTailList=[]
        if aTail in sequence:
            taillength=0
            reverseSequence = str(Seq(sequence).reverse_complement())
            findings = [int(x) for x in find_sequence(reverseSequence, tTail)]
            df2.at[row,'Tail_Seed_Hits']=len(findings)
            if len(findings)>0:
                aTailFlag=1
                amin = min(findings)
                iterable = findings
                groupings = [list(group) for group in mit.consecutive_groups(iterable)]
                flag=0
                for group in groupings:
                    if flag==0:

                        taillength+= len(aTail) + (len(group)-1)
                        groupEnd = max(group)+len(aTail)
                        flag=1

                    else:

                        if min(group)-groupEnd<=4:
                            taillength+= (len(aTail) + (len(group)-1) + abs(min(group)-groupEnd))
                            groupEnd = max(group)+len(aTail)
                        else:
                            continue

                aTailLength = taillength
                
            else:
                pass
        
        if tTailFlag>0 and aTailFlag==0:
            
            df2.at[row,'Tail_Begins']=tmin
            df2.at[row,'Tail_Type']='Possible_T-Tail'
            df2.at[row,'Tail_Length']=tTailLength
            
        elif tTailFlag==0 and aTailFlag>0:
            
            df2.at[row,'Tail_Begins']=amin
            df2.at[row,'Tail_Type']='Possible_A-Tail'
            df2.at[row,'Tail_Length']=aTailLength
        
        
        elif aTailFlag>0 and tTailFlag>0:
            
            orientation = str(df2.at[row,'Orientation'])
            
            if tTailLength<aTailLength and amin<tmin:
                df2.at[row,'Tail_Begins']=amin
                df2.at[row,'Tail_Type']='Possible_A-Tail*_and_Possible_T-Tail'
                df2.at[row,'Tail_Length']=aTailLength
                
            elif tTailLength>aTailLength and tmin<amin:
                df2.at[row,'Tail_Begins']=tmin
                df2.at[row,'Tail_Type']='Possible_A-Tail_and_Possible_T-Tail*'
                df2.at[row,'Tail_Length']=tTailLength
                
            else:
                if orientation !='None':
                    
                    if tmin>=50 and amin<50:
                        df2.at[row,'Tail_Begins']=amin
                        df2.at[row,'Tail_Type']='Possible_A-Tail*_and_Possible_T-Tail'
                        df2.at[row,'Tail_Length']=aTailLength
                    elif tmin<50 and amin>=50:
                        df2.at[row,'Tail_Begins']=tmin
                        df2.at[row,'Tail_Type']='Possible_A-Tail_and_Possible_T-Tail*'
                        df2.at[row,'Tail_Length']=tTailLength
                    
                        
                    elif aTailLength<tTailLength and orientation == '-':
                        df2.at[row,'Tail_Begins']=tmin
                        df2.at[row,'Tail_Type']='Possible_A-Tail_and_Possible_T-Tail*'
                        df2.at[row,'Tail_Length']=tTailLength

                    elif aTailLength>tTailLength and orientation == '+':
                        df2.at[row,'Tail_Begins']=amin
                        df2.at[row,'Tail_Type']='Possible_A-Tail*_and_Possible_T-Tail'
                        df2.at[row,'Tail_Length']=aTailLength
                        
                    elif orientation == '+':
                        df2.at[row,'Tail_Begins']=amin
                        df2.at[row,'Tail_Type']='Possible_A-Tail*_and_Possible_T-Tail'
                        df2.at[row,'Tail_Length']=aTailLength
                    
                    
                    elif orientation == '-':
                        df2.at[row,'Tail_Begins']=tmin
                        df2.at[row,'Tail_Type']='Possible_A-Tail_and_Possible_T-Tail*'
                        df2.at[row,'Tail_Length']=tTailLength
                    
                    else:
                        continue
                else:
                    continue
                    
                
        
        else:
            continue
            
    return(df2)
            
   
def aluLinker(df):
    df2 = df.copy()
    
    for row in df2.index:
        
        threePrimeFlags=0
        
        if df2.at[row,'TE_Designation'] == 'SINE/Alu' and str(df2.at[row,'Tail_Type'])!='No_Tail_Type':
            
            tailDesignation = str(df2.at[row,'Tail_Type'])
            tailStartSite = int(df2.at[row,'Tail_Begins'])
            elements = ast.literal_eval(str(df2.at[row,'TE_Hits']))
            
            if len(elements)==1:
                for element in elements:

                    splitElement = element.split()

                    if abs(int(splitElement[12])-tailStartSite)>=120 and abs(int(splitElement[12])-tailStartSite)<=150:
                        threePrimeFlags+=1
                    else:
                        continue
            else:
                pass
                        
            if threePrimeFlags>0:
                df2.at[row,'FILTER_RESULTS']=['Alu_Linker_Region_Warning']
            else:
                continue
            
            
        else:
            continue
    return(df2) 



def finalQuickCheck(df):
    
    df2 = df.copy()
    
    for row in df2.index:
            
            if str(df2.at[row,'TE_Designation']) =='LINE/L1':
                if str(df2.at[row,'Element_Annotation']) in [youngElement for youngElement in str(args.youngLINE).split(",")]:
                    continue
                else:
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['OLDER_LINE_SUBFAMILY']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('OLDER_LINE_SUBFAMILY')
                        df2.at[row,'FILTER_RESULTS']=newList
                    
                    
            elif str(df2.at[row,'TE_Designation']) =='SINE/Alu':
                
                aluList = ast.literal_eval(str(df2.at[row,'TE_Hits']))
                if len(aluList)==1:
                    continue
                else:
                    
                    aluHitList=[]
                    for hit in aluList:
                        splitHit = hit.split()
                        if 'Alu' in str(splitHit[9]):
                            aluHitList.append(splitHit[9])
                        else:
                            continue
                            
                    if len(set(aluHitList))>1:
                        if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                            df2.at[row,'FILTER_RESULTS']=['MULTI-ALU_Hits']
                        else:
                            newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                            newList.append('MULTI-ALU_Hits')
                            df2.at[row,'FILTER_RESULTS']=newList

                    else:
                        continue

    return(df2)


def finalQuickCheck2(df):
    
    df2 = df.copy()
    
    for row in df2.index:
            
            if str(df2.at[row,'TE_Designation']) =='LINE/L1':
                if int(df2.at[row,'Sequence_Length']) >int(str(args.Lengths).split(",")[1]):
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['LINE_>Len']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('LINE_>10kLen')
                        df2.at[row,'FILTER_RESULTS']=newList
                else:
                    pass
                    
                if int(df2.at[row,'Element_Divergence']) >float(str(args.Divergences).split(",")[1]):
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['LINE_DIVERGENCE']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('LINE_DIVERGENCE')
                        df2.at[row,'FILTER_RESULTS']=newList
                else:
                    continue
                    
                    
            elif str(df2.at[row,'TE_Designation']) =='SINE/Alu':
                
                if int(df2.at[row,'Sequence_Length']) > int(str(args.Lengths).split(",")[0]):
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['ALU_>Len']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('ALU_>500Len')
                        df2.at[row,'FILTER_RESULTS']=newList
                else:
                    pass
                    
                if int(df2.at[row,'Element_Divergence']) >float(str(args.Divergences).split(",")[0]):
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['ALU_DIVERGENCE']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('ALU_DIVERGENCE')
                        df2.at[row,'FILTER_RESULTS']=newList
                else:
                    continue
                        
                        
            elif str(df2.at[row,'TE_Designation']) =='Retroposon/SVA':
                if int(df2.at[row,'Sequence_Length']) >int(str(args.Lengths).split(",")[2]):
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['SVA_>Len']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('SVA_>10kLen')
                        df2.at[row,'FILTER_RESULTS']=newList
                else:
                    pass
                    
                if int(df2.at[row,'Element_Divergence']) >float(str(args.Divergences).split(",")[2]):
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['SVA_DIVERGENCE']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('SVA_DIVERGENCE')
                        df2.at[row,'FILTER_RESULTS']=newList
                else:
                    continue
            
            else:
                continue

    return(df2)


def tailCounterCheck(df):
    df2 = df.copy()
    
    for row in df2.index:
        if df2.at[row,'Tail_Type']!='No_Tail_Type':
            flag=0

            if int(df2.at[row,'Tail_Begins'])/int(df2.at[row,'Sequence_Length'])>.5 or int(df2.at[row,'Tail_Begins'])>int(args.maxTailStart):
                flag+=1
            else:
                pass

            if flag ==0:
                continue
            else:
                flagg='Bad_Tail_Position'
                if df2.at[row,'FILTER_RESULTS']!='Good_Row':
                    tempFlagList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                    tempFlagList.append(flagg)
                    df2.at[row,'FILTER_RESULTS'] = tempFlagList
                else:
                    df2.at[row,'FILTER_RESULTS'] = [flagg]
        else:
            continue
    
    return(df2)



def cleanDataTEPercentage(df):
    
    df2 = df.copy() 
    df2['Element_Annotation']='No_Element_Annotation'
    df2['Element_Divergence']=0.0
    df2['TE_Proportion']='NONE'
    df2['TE_Percentage']=0.0
    
    #For each locus
    for row in tqdm(df2.index):
        
        if df2.at[row,'TE_Hits'] == 'NONE' or int(df2.at[row,'Sequence_Length'])>int(args.Upper):
            continue
        else:
            

            columnName = 'Sequence'

            tempDict={int(x):0 for x in range(1,len(df2.at[row,columnName])+1)}    
            tempDict2={}
            tempDict3={}
            tempDivDict={}

            teHitList = ast.literal_eval(str(df2.at[row,'TE_Hits']))

            for hit in teHitList:

                splitHit = hit.split()
                focusElement = str(splitHit[9])
                tempDict3[focusElement]=str(splitHit[10])

                if focusElement in tempDict2.keys():
                    tempDivDict[focusElement].append(float(splitHit[1]))
                    pass
                else:
                    tempDivDict[focusElement]=[float(splitHit[1])]
                    tempDict2[focusElement]={x:0 for x in range(1,len(df2.at[row,columnName])+1)}
                    pass

                for coordinate in range(int(splitHit[5]), int(splitHit[6])+1):
                    tempDict[coordinate]+=1
                    tempDict2[focusElement][coordinate]+=1

            tePercentage = len([x for x in tempDict.values() if x >0])/len(tempDict)

            df2.at[row,'TE_Percentage']=float(tePercentage)

            #print(tempDict3)
            #print(tempDivDict)

            if len(tempDict2)==1:
                mykey = [x for x in tempDict2.keys()][0]
                df2.at[row,'TE_Designation']=str(tempDict3[mykey])
                df2.at[row,'TE_Proportion']={mykey:tePercentage}
                df2.at[row,'Element_Annotation']=mykey
                df2.at[row,'Element_Divergence']=np.median(tempDivDict[mykey][0])

            else:

                tempDict4 = {x:len([y for y in tempDict2[x].values() if y>0])/len(tempDict) for x in tempDict2.keys()}
                #print(tempDict4)
                maxKey = max(tempDict4, key=tempDict4.get)
                df2.at[row,'TE_Designation']=tempDict3[str(maxKey)]
                df2.at[row,'TE_Proportion']=tempDict4
                df2.at[row,'Element_Annotation']=maxKey
                df2.at[row,'Element_Divergence']=np.median(tempDivDict[maxKey])
        
    return(df2)


def repeatmaskerPatternFilter(df):
    df2 = df.copy()
    annotationList=[]
    df2['Unique_Element_Count']='One_Element'
    for row in df2.index:
        if df2.at[row,'TE_Proportion'] == 'NONE':
            continue
        else:

            prodict = ast.literal_eval(str(df2.at[row,'TE_Proportion']))

            if str(df2.at[row,'TE_Designation']) == 'SINE/Alu':
                aluList = [x for x in prodict.keys() if 'ALU' in str(x).upper()]
                if len(aluList)==1:
                    continue
                else:
                    df2.at[row,'Unique_Element_Count']="More_Than_One_Element"
            else:
                continue 
    
    
    for row in df2.index:
        if df2.at[row,'TE_Proportion'] == 'NONE':
            continue
        else:
        
            if df2.at[row,'Unique_Element_Count'] == 'One_Element':

                tempDict = ast.literal_eval(str(df2.at[row,'TE_Proportion']))
                dictLength = len([x for x in tempDict.keys()])

                if str(df2.at[row,'TE_Designation']) == 'SINE/Alu':
                    if len(ast.literal_eval(str(df2.at[row,'TE_Hits'])))>dictLength:
                        df2.at[row,'Unique_Element_Count'] = 'One_Element_ODD'
                    else:
                        continue



                elif str(df2.at[row,'TE_Designation']) == 'LINE/L1':
                    if len(ast.literal_eval(str(df2.at[row,'TE_Hits'])))>dictLength:
                        df2.at[row,'Unique_Element_Count'] = 'One_Element_ODD'
                    else:
                        continue

                else:
                    if len(ast.literal_eval(str(df2.at[row,'TE_Hits'])))>dictLength:
                        df2.at[row,'Unique_Element_Count'] = 'One_Element_ODD'
                    else:
                        continue


            else:
                continue
    return(df2)


def findTwinPriming(df):
    df2 = df.copy()
    df2['Twin_Priming_Flag']='NONE'
    
    for row in df2.index:
        if df2.at[row,'TE_Designation']=='LINE/L1':
            allhits = ast.literal_eval(str(df2.at[row,'TE_Hits']))
            orientations = []
            for hit in allhits:
                splitHit=hit.split()
                if str(splitHit[10]) == 'LINE/L1':
                    orientations.append(splitHit[8])
                else:
                    continue
            if len(set(orientations))>1:
                df2.at[row,'Twin_Priming_Flag']='FLAG'
            else:
                continue
            
            
        else:
            continue
    return(df2)


def simpleRepeatCheck(df):
    df2 = df.copy()
    
    for row in df2.index:
        if df2.at[row,'TE_Designation']=='Simple_repeat':
            tempDict = {x:float(y) for x,y in ast.literal_eval(str(df2.at[row,'TE_Proportion'])).items() if ')n' not in x}
            if len(tempDict)>0:
                maxKey = str(max(tempDict, key=tempDict.get)).upper()
                if 'ALU' in maxKey and len(tempDict) == 1:
                    #print(row)
                    df2.at[row,'TE_Designation']= 'SINE/Alu'
                    df2.at[row,'Element_Annotation']= maxKey

                elif'L1' in maxKey and len(tempDict) == 1:
                    #print(row)
                    df2.at[row,'TE_Designation']= 'LINE/L1'
                    df2.at[row,'Element_Annotation']= maxKey

                elif 'SVA' in maxKey and len(tempDict) == 1:
                    #print(row)
                    df2.at[row,'TE_Designation']= 'Retroposon/SVA'
                    df2.at[row,'Element_Annotation']= maxKey

                else:
                    continue
            else:
                continue
        else:
            continue
    return(df2)

def main():

	global args 
	# Ask for user inputs
	parser = ap.ArgumentParser()
	parser.add_argument("-i", "--input", dest="input", required=True,
	                help="Please specify the path to the fasta file that you ran RepeatMasker on.")

	parser.add_argument("-r", "--repeats", dest="Repeat", required=True,
	                help="Please specify the path to the RepeatMasker output file.")

	parser.add_argument("-rd", "--RepeatDiv", dest='RepeatDiv', default=20)

	parser.add_argument("-d", "--divergences", dest='Divergences', default='6,15,15')

	parser.add_argument("-l", "--lengths", dest='Lengths', default='500,10000,10000')

	parser.add_argument("-u", "--upper", dest='Upper', default=50000)

	parser.add_argument("-t", "--maxTailStart", dest='maxTailStart', default=50)

	parser.add_argument("-yl", "--youngLINE", dest='youngLINE', default='L1HS,L1PA2')

	parser.add_argument("-o", "--output", dest='Output', required=True, help="Please specify the output path and filename for the csv file")

	args = parser.parse_args()

	sequences ={}

	fasta_sequences = SeqIO.parse(open(str(args.input)),'fasta')
	for fasta in fasta_sequences:
	    name, sequence = fasta.id, str(fasta.seq)
	    sequences[name]=sequence
	insDF = pd.DataFrame(data=[[x,y] for x,y in sequences.items()], columns=['ID','Sequence']).set_index("ID")

	rmdirectory = str(args.Repeat)
	callDict={}
	with open(rmdirectory, 'r') as file:
	    lines_after_header = file.readlines()[3:]
	    for line in lines_after_header:
	        goodline = ' '.join(line.split())
	        if str(goodline.split(" ")[4]) in callDict.keys() and float(goodline.split(" ")[1])<=float(args.RepeatDiv):
	            callDict[str(goodline.split(" ")[4])]['Annotations'].append(goodline)
	        elif float(goodline.split(" ")[1])<=float(args.RepeatDiv):
	            callDict[str(goodline.split(" ")[4])]={'Annotations':[]}
	            callDict[str(goodline.split(" ")[4])]['Annotations'].append(goodline)
	        else:
	            continue
	file.close()
	insDF['TE_Hits']='NONE'
	for row in insDF.index:
	    if row in callDict.keys():
	        insDF.at[row,'TE_Hits']=callDict[row]['Annotations']
	    else:
	        continue

	insDF['Sequence_Length']=[len(insDF.at[x,'Sequence']) for x in insDF.index]


	#Run the functions (I know this is ugly, but it makes it simpler to debug when the functions aren't nested)
	print('Reading in the RepeatMasker info')
	insDF_Filtered= orientationFinder(cleanDataTEPercentage(insDF))
	print('Looking for tails')
	insDF_Filtered2 = tailCounter(insDF_Filtered)
	print('Checking if the tail is actually an Alu A-linker region')
	insDF_Filtered3 = aluLinker(insDF_Filtered2)
	print('Checking tail results')
	insDF_Filtered4 = tailCounterCheck(insDF_Filtered3)
	print('Running some length, divergence checks')
	insDF_Filtered5 =finalQuickCheck(insDF_Filtered4)
	insDF_Filtered6 = finalQuickCheck2(insDF_Filtered5)
	insDF_Filtered7 = repeatmaskerPatternFilter(insDF_Filtered6)
	print('Checking for potential signs of twin priming for L1s')
	insDF_Filtered8 = findTwinPriming(insDF_Filtered7)
	print('Final checks.')
	insDF_Filtered9 = simpleRepeatCheck(insDF_Filtered8)

	print('Writing out the results. Remember to double-check things yourself as this program is in beta!')
	insDF_Filtered9.to_csv(str(args.Output))



if __name__ == "__main__":
        main()