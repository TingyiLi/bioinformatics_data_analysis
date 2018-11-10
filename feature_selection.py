#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 13:57:48 2018

@author: litingyi
"""
import pandas as pd
import numpy as np
import sys
import seaborn as sns
import os
import re
import math
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import Normalizer
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import SelectKBest
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from minepy import MINE
from numpy import random
from sklearn.linear_model import LogisticRegression
from sklearn import metrics

def read_data(pathway):
    all_features=pd.DataFrame([])
    for dir in os.listdir(pathway):
        try:
            result=pd.read_csv(pathway+dir)
            id_num=int(re.sub('\D','',dir))
            result.insert(0,'Patient ID',id_num)
            all_features=pd.concat([all_features,result],axis=0)
        except:
            print (dir,"is empty!")
    return all_features

#和求相关系数方法类似,由于MINE的设计不是函数式的，定义mic方法将其为函数式的，返回一个二元组，二元组的第2项设置成固定的P值0.5
#pic=1-mic
def pic(x,y):
    m = MINE()
    m.compute_score(x,y)
    return [1-m.mic(),0.5]

def plot_corr_p(corr_pvalue):
    plt.figure(figsize=(7,5))
    plt.subplot(211)
    plt.plot(corr_pvalue[0])
    plt.title('correlation')
    plt.subplot(212)
    plt.plot(corr_pvalue[1])
    plt.yscale('log')
    plt.title('p_value')
    return

def best_k_features(corr,k):
    corr[0]=map(lambda x: abs(float(x)),corr[0])
    feature_corr=np.array(corr).T.tolist()
    feature_corr.sort(key=lambda x:x[0],reverse=True)
    return np.array(feature_corr).T[2][0:k-1].tolist(),np.array(feature_corr).T[0][0:k-1].tolist()

"""build imblance-adjusted bootstrap training set and test set"""
def buildBootSet(Y,nBoot):
    idPatients=Y.index.values.tolist()
    idPos=Y[Y['LungMets']==1].index.values  # patient id of positive instances
    idNeg=Y[Y['LungMets']==0].index.values  # patient id of negative instances
    nPos=len(idPos)
    nNeg=len(idNeg)
    nInst=nPos+nNeg
    idPos=np.tile(idPos,(nNeg,1)).flatten().tolist()
    idNeg=np.tile(idNeg,(nPos,1)).flatten().tolist()
    idInst=idPos+idNeg
    random.shuffle(idInst)
    random.shuffle(idInst)
    nSize=len(idInst)
    randIndex=np.random.randint(0,nSize,size=(nInst,nBoot))
    bootsam=[np.array(idInst)[randIndex[:,i]] for i in range(nBoot)]
    testSets=[list(set(idPatients)-set(bootsam[i])) for i in range(nBoot)]
    """verification for each bootstrap smaple, both training and testing sets must have at least 2 instances of each class, and these instances must be different"""
    for n in range(nBoot):
        training=Y.loc[bootsam[n],:]
        testing=Y.loc[testSets[n],:]
        while (len(set(training[training['LungMets']==0].index.values.tolist()))<2 or len(set(training[training['LungMets']==1].index.values.tolist()))<2 or len(set(testing[testing['LungMets']==1].index.values.tolist()))<2 or len(set(testing[testing['LungMets']==0].index.values.tolist()))<2):
            bootsam[n]=np.array(idInst)[np.random.randint(0,nSize,nInst)]
            testSets[n]=list(set(idPatients)-set(bootsam[n]))
            training=Y.loc[bootsam[n],:]
            testing=Y.loc[testSets[n],:]
    return bootsam,testSets


def Cal_spearman_bootstrap(X,Y,bootSam):
    B=len(bootSam) #bootstrap 次数
    matData=np.array([X.loc[bootSam[b],:].values for b in range(B)]) #所有bootstrap产生的样本，大小为B*50*nFeatures
    nFeatures=matData.shape[2]
    matCorr=np.zeros(nFeatures)

    for b in range(B):
        dataBoot=matData[b]
        outcomeBoot=np.array(Y.loc[bootSam[b],:].values)
        spearman_cal=lambda X,Y: np.array(map(lambda x: spearmanr(x,Y), X.T)) #计算spearman系数
        spearman_corr=spearman_cal(dataBoot,outcomeBoot)
        matCorr=matCorr+spearman_corr[:,0]

    matCorr=abs(matCorr/B)

    return matCorr

def Cal_PIC_bootstrap(X,Yname,bootSam):
    B=len(bootSam) #bootstrap 次数
    XData=np.array([X.loc[bootSam[b],:].values for b in range(B)]) #所有bootstrap产生的样本，大小为B*50*otherFeaturesNum
    FeaturesNum=XData.shape[2]
    Y=pd.DataFrame(X.loc[:,Yname])
    YData=np.array([Y.loc[bootSam[b],:].values for b in range(B)]) #所有bootstrap产生的样本，大小为B*50*1
    matPIC=np.zeros(FeaturesNum)
    for b in range(B):
        XBoot=XData[b]
        YBoot=YData[b][:,0]
        pic_score=lambda X,Y:np.array(map(lambda x: pic(x,Y),X.T)) #计算pic系数
        pic_corr=pic_score(XBoot,YBoot)
        matPIC=matPIC+pic_corr[:,0]
    matPIC=matPIC/B

    return matPIC

def featureReduction(featureData,outcomeData,nBoot,alpha,delta,setSize):
    fSet=[] #featureReduction 函数返回值，返回feature的name
    """Getting imblance-adjusted bootstrap samples for Gain equation"""
    bootsam,testSets=buildBootSet(outcomeData,nBoot)
    #featureData=pd.DataFrame(MinMaxScaler().fit_transform(featureData),columns=featureData.columns)  #无量纲化--归一化（二范数）

    """Geting bootstrap results for the first part of the Gain equation"""
    matCorr=Cal_spearman_bootstrap(featureData,outcomeData,bootsam)

    """Choosing first feature (depends only on spearman correlation)"""
    Corr_df=pd.DataFrame(matCorr,index=featureData.columns,columns=['spearman_corr'])
    first_feature=Corr_df[Corr_df['spearman_corr']==max(Corr_df['spearman_corr'])].index[0]
    print (first_feature,max(Corr_df['spearman_corr']))
    fSet.append(first_feature)
    Corr_df=Corr_df.drop(first_feature)

    """Computing for other features"""
    allFeaturesNum=featureData.values.shape[1]
    PICtest=np.zeros((allFeaturesNum,setSize-1))
    for f in range(1,setSize):
        selectedFeatureName=fSet[f-1]
        PICtest[:,f-1]=PICtest[:,f-1]+Cal_PIC_bootstrap(featureData,selectedFeatureName,bootsam)
        PICtemp = np.zeros(allFeaturesNum)
        for k in range(1,f):
            PICtemp=PICtemp+2*(f-k+1)/(f*(f+1))*PICtest[:,k-1]
        matPic=pd.DataFrame(PICtemp,index=featureData.columns,columns=['PIC']).drop(fSet).values[:,0]
        matCorr=Corr_df.values[:,0]
        Gain=alpha*matCorr+delta*matPic
        indChosen=np.argwhere(Gain == max(Gain))[0]
        featureChosen=Corr_df.index[indChosen][0]
        print (featureChosen,max(Gain))
        fSet.append(featureChosen)
        Corr_df=Corr_df.drop(featureChosen)

    return fSet

def featureSelection(X,Y,maxOrder,nBoot):
    models=[] #输出
    models_auc=[]
    nFeat=X.values.shape[1]
    """Getting imblance-adjusted bootstrap samples for Logistic regression"""
    trainSets,testSets=buildBootSet(Y,nBoot)

    """Forward feature selection  (for all starters)"""
    modelMat=[['' for i in range(maxOrder)] for j in range(nFeat)]
    metricMat=np.zeros((nFeat,maxOrder))
    for i in range(nFeat):
        indLeft=X.columns.values.tolist()
        #order 1
        starter=[indLeft[i]]
        print (starter,i)
        modelMat[i][0]=indLeft[i]
        del indLeft[i]
        aucBoot,sensBoot,sepcBoot=bootstrap_AUC_sens_spec(trainSets,testSets,X,Y,starter)
        metricMat[i,0]=1*aucBoot+0*(1-abs(sensBoot-sepcBoot))

        #Going for orders 2 to maxOrder
        for j in range(1,maxOrder):
            maxMetric=0
            for k in range(nFeat-j):
                indexModel=modelMat[i][0:j]+[indLeft[k]]
                aucBoot,sensBoot,sepcBoot=bootstrap_AUC_sens_spec(trainSets,testSets,X,Y,indexModel)
                metricTemp=1*aucBoot+0*(1-abs(sensBoot-sepcBoot))
                if metricTemp>=maxMetric:
                    maxMetric=metricTemp
                    index=indLeft[k]

            modelMat[i][j]=index
            metricMat[i,j]=maxMetric
            indLeft.remove(index)

    #Obtaining maximum AUC for every model order
    for i in range(maxOrder):
        current_index=np.argwhere(metricMat[:,i]==max(metricMat[:,i]))[0][0]
        featureName=modelMat[current_index][0:i+1]
        models.append(featureName)
        models_auc.append(max(metricMat[:,i]))

    return models,models_auc

def bootstrap_AUC_sens_spec(trainSets,testSets,X,Y,features):
    top=1-1/np.exp(1)
    low=1/np.exp(1)
    nBoot=len(trainSets)
    Xtrain=X.loc[:,features]
    Xtrain=MinMaxScaler().fit_transform(Xtrain)
    clf_lr=LogisticRegression()
    clf_lr.fit(Xtrain,Y)
    predic_prob_y=clf_lr.predict_proba(Xtrain)
    labels=clf_lr.predict(Xtrain)
    sensData,specData=cal_sens_spec(labels,Y)
    aucData=metrics.roc_auc_score(Y,predic_prob_y[:,1])
    aucTemp=0;sensTemp=0;specTemp=0
    for n in range(nBoot):
        Xtrain=X.loc[trainSets[n],features];Xtest=X.loc[testSets[n],features];Ytrain=Y.loc[trainSets[n],:];Ytest=Y.loc[testSets[n],:]
        Xtrain=MinMaxScaler().fit_transform(Xtrain);Xtest=MinMaxScaler().fit_transform(Xtest)
        clf_lr=LogisticRegression()
        clf_lr.fit(Xtrain,Ytrain)
        predic_prob_y=clf_lr.predict_proba(Xtest)
        labels=clf_lr.predict(Xtest)
        sensBoot,specBoot=cal_sens_spec(labels,Ytest)
        aucBoot=metrics.roc_auc_score(Ytest,predic_prob_y[:,1])

        #For AUC
        alpha=top/(1-low*(aucData-aucBoot)/(aucData-0.5+np.spacing(1)))
        if alpha>1:
            alpha=1
        elif alpha<top:
            alpha=top
        if aucBoot<0.5:
            aucBoot=0.5
        aucTemp=aucTemp+(1-alpha)*aucData+alpha*aucBoot

        #For sensitivity
        alpha=top/(1-low*(sensData-sensBoot)/sensData+np.spacing(1))
        if alpha<top:
            alpha=top
        sensTemp=sensTemp+(1-alpha)*sensData+alpha*sensBoot

        #For specificity
        alpha=top/(1-low*(specData-specBoot)/specData+np.spacing(1))
        if alpha<top:
            alpha=top
        specTemp=specTemp+(1-alpha)*specData+alpha*specBoot
    aucTemp=aucTemp/nBoot;sensTemp=sensTemp/nBoot;specTemp=specTemp/nBoot
    return aucTemp,sensTemp,specTemp

def cal_sens_spec(predLabel,actuLabel):
    Label=pd.DataFrame(actuLabel.values,columns=['actual_label'])
    Label['predicted_label']=predLabel
    tp=len(Label[(Label['actual_label']==1)&(Label['predicted_label']==1)])
    tn=len(Label[(Label['actual_label']==0)&(Label['predicted_label']==0)])
    fp=len(Label[(Label['actual_label']==0)&(Label['predicted_label']==1)])
    fn=len(Label[(Label['actual_label']==1)&(Label['predicted_label']==0)])
    sensitivity=float(tp)/(tp+fn)
    specificity=float(tn)/(tn+fp)
    return sensitivity,specificity


def UnivariateAnalysis(X,Y):
    #X=pd.DataFrame(MinMaxScaler().fit_transform(X),columns=data.columns[0:-1])  #无量纲化--归一化（二范数）
    features=np.array(X.values)
    #VarianceThreshold(threshold=3).fit_transform(features)
    #cal_corr=lambda X,Y: np.array(map(lambda x: pearsonr(x,Y), X.T)) #计算pearson系数，这里会报错，是不是适用？？？？
    cal_corr=lambda X,Y: np.array(map(lambda x: spearmanr(x,Y), X.T)) #计算spearman系数
    #mine_score=lambda X,Y:np.array(map(lambda x: mic(x,Y),X.T))  #计算互信息系数,有UnicodeDecodeError:
    label=np.array(data.loc[:,['LungMets']])
    corr_pvalue=cal_corr(features,label).T.tolist()
    plot_corr_p(corr_pvalue)
    #corr_pvalue.append(data.columns.values.tolist()[0:-1])
    #best_25_features,best_25_corr=best_k_features(corr_pvalue,25)
    features_Heatmap(X,Y)
    return

def features_Heatmap(X,Y):
    """画出feature_reduction 和输出的heatmap"""
    data_reduction=X
    data_reduction['LungMets']=Y.loc[:,['LungMets']]
    plt.figure(figsize=(12,10))
    corr = sns.heatmap(data_reduction.corr(),vmax=0.8,annot=True)
    return

def evaluate_model(models,X,Y,nBoot):
    trainSets,testSets=buildBootSet(Y,nBoot)
    orders=len(models)
    evalutes_df=pd.DataFrame([],columns=['values','kind','Order'])
    for r in range(orders):
        for n in range(nBoot):
            Xtrain=X.loc[trainSets[n],models[r]];Xtest=X.loc[testSets[n],models[r]];Ytrain=Y.loc[trainSets[n],:];Ytest=Y.loc[testSets[n],:]
            Xtrain=MinMaxScaler().fit_transform(Xtrain);Xtest=MinMaxScaler().fit_transform(Xtest)
            clf_lr=LogisticRegression()
            clf_lr.fit(Xtrain,Ytrain)
            predic_prob_y=clf_lr.predict_proba(Xtest)
            labels=clf_lr.predict(Xtest)
            sensBoot,specBoot=cal_sens_spec(labels,Ytest)
            aucBoot=metrics.roc_auc_score(Ytest,predic_prob_y[:,1])
            temp_df=pd.DataFrame([[aucBoot,'AUC',r+1],[sensBoot,'Sensitivity',r+1],[specBoot,'Specitivity',r+1]],columns=['values','kind','Order'])
            evalutes_df=pd.concat((evalutes_df,temp_df),axis=0)
    sns.set(style="whitegrid")
    g = sns.factorplot(x="Order", y="values", hue="kind", data=evalutes_df,capsize=.2, palette="YlGnBu_d", size=6, aspect=.75)

    return

if __name__=='__main__':
    """读取数据文件"""
    pathway='./features/CT/'
    data=read_data(pathway)
    outcome_ori=pd.read_csv('./outcome_vector.csv')
    data=pd.merge(data,outcome_ori,on='Patient ID',how='inner')

    """分离特征和输出"""
    data.drop(['general_info_BoundingBox','general_info_EnabledImageTypes','general_info_GeneralSettings','general_info_ImageHash','general_info_ImageSpacing','general_info_MaskHash','general_info_Version','general_info_VolumeNum'],axis=1,inplace=True)
    data.index=data['Patient ID']
    data=data.drop(['Patient ID'],axis=1)
    X=pd.DataFrame(data.iloc[:,0:-1] )#自变量
    Y=pd.DataFrame(data.iloc[:,[-1]])  #因变量

    """
    #依据STS文献计算bootstrap Gain，进而筛选25个特征
    features_25=featureReduction(X,Y,500,0.5,0.5,25)
    #np.save('./bootstrap_reduction_25features_1.npy',features_25)
    X_reduction=pd.DataFrame(X.loc[:,features_25])
    UnivariateAnalysis(X_reduction,Y)
    sys.exit()

    #一元分析filter特征，依据方差和相关性
    UnivariateAnalysis(X,Y)
    """

    #利用sklearn构建LR模型，求解AUC
    features=np.load('./bootstrap_reduction_25features_1.npy')
    features=[x.decode('utf-8') for x in features]
    X_reduction=pd.DataFrame(X.loc[:,features])
    modelsByOrders,models_auc=featureSelection(X_reduction,Y,10,1000)
    np.save('./modelsByOrder_200bootstrap.npy',modelsByOrders)
    np.save('./models_auc_200bootstrap.npy',models_auc)
    evaluate_model(modelsByOrders,X_reduction,Y,10)
