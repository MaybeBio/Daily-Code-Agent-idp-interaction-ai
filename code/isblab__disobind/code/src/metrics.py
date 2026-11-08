########## Calculate evaluation metrics using torchmetrics #########
########## ------>"May the Force serve u well..." <------###########
####################################################################

############# One above all #############
##-------------------------------------##np
import torch
from torchmetrics.classification import ( BinaryStatScores, BinaryAUROC, BinaryAveragePrecision,
                                            BinaryRecall, BinaryPrecision, BinaryF1Score, BinaryAccuracy,
                                            BinaryMatthewsCorrCoef )


## Calculate recall, precision, f1-score, accuracy, mcc ##
##########################################################
def torch_metrics( preds, target, threshold, multidim_avg, device ):
    """
    Calculates the following evaluation metrics:
        Recall gives a measure of how many true positives have been predicted
          out of the actual no. of positives
        Precision gives  a measure of hoe many true positives have been predicted
            out of the total no. of predicted positives
        F1-score is  measure of a models accuracy
            It's the harmonic mean of precision and recall
        Accuracy is a measure of how far or close are the predicted results from the actual values
            Misleading for imbalanced datasets
    """
    target = torch.where( target > threshold, 1, 0 )
    # multidim_average = "samplewise" # samplewise, global
    m_avg = multidim_avg.split( "_" )[0]

    recall    = BinaryRecall(    threshold = threshold, multidim_average = m_avg ).to( device )
    precision = BinaryPrecision( threshold = threshold, multidim_average = m_avg ).to( device )
    f1        = BinaryF1Score(   threshold = threshold, multidim_average = m_avg ).to( device )
    mcc       = BinaryMatthewsCorrCoef(  threshold = threshold ).to( device )
    bap       = BinaryAveragePrecision( thresholds = 10 ).to( device )
    acc       = BinaryAccuracy(  threshold = threshold, multidim_average = m_avg ).to( device )
    auroc     = BinaryAUROC( thresholds = 10 ).to( device )

    if multidim_avg == "global":
        return [recall( preds, target ), precision( preds, target ), f1( preds, target ), bap( preds, target ), 
                mcc( preds, target ), auroc( preds, target ), acc( preds, target )]

    elif multidim_avg == "samplewise":
        recall = torch.mean( recall( preds, target ) )
        precision = torch.mean( precision( preds, target ) )
        f1 = torch.mean( f1( preds, target ) )
        bap = torch.mean( bap( preds, target ) )
        mcc = torch.mean( mcc( preds, target ) )
        auroc = torch.mean( auroc( preds, target ) )
        acc = torch.mean( acc( preds, target ) )
        return [recall, precision, f1, bap, mcc, auroc, acc]

    elif multidim_avg == "samplewise_none":
        recall = recall( preds, target )
        precision = precision( preds, target )
        f1 = f1( preds, target )
        bap = bap( preds, target )
        mcc = mcc( preds, target )
        auroc = auroc( preds, target )
        acc = acc( preds, target )
        return [recall, precision, f1, bap, mcc, auroc, acc]
