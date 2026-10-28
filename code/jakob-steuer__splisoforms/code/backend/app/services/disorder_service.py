import numpy as np
from sqlalchemy.orm import Session
from sqlalchemy import text
from aiupred import AIUPred
from app.models.bio import DisorderScore, Isoform
import logging

logger = logging.getLogger(__name__)

class DisorderService:
    _predictor = None

    @classmethod
    def get_predictor(cls):
        if cls._predictor is None:
            logger.info("Initializing AIUPred predictor (CPU-mode fallback)...")
            cls._predictor = AIUPred()
        return cls._predictor

    def analyze_isoform(self, db: Session, isoform_id: int) -> DisorderScore:
        # 1. Check if already exists
        existing = db.query(DisorderScore).filter(DisorderScore.isoform_id == isoform_id).first()
        if existing:
            return existing

        # 2. Get sequence
        isoform = db.query(Isoform).filter(Isoform.id == isoform_id).first()
        if not isoform or not isoform.sequence:
            raise ValueError(f"Isoform {isoform_id} not found or has no sequence.")

        sequence = isoform.sequence
        predictor = self.get_predictor()

        # 3. Run predictions
        try:
            disorder_values = predictor.predict_disorder(sequence)
            binding_values = predictor.predict_binding(sequence)
            
            # Convert numpy arrays to lists for SQLAlchemy ARRAY support
            disorder_list = disorder_values.tolist()
            binding_list = binding_values.tolist()

            new_score = DisorderScore(
                isoform_id=isoform_id,
                disorder_scores=disorder_list,
                anchor_scores=binding_list
            )
            db.add(new_score)
            db.commit()
            db.refresh(new_score)
            return new_score
        except Exception as e:
            db.rollback()
            logger.error(f"Error predicting disorder for isoform {isoform_id}: {str(e)}")
            raise e

    def get_disorder_data(self, db: Session, isoform_id: int):
        score = db.query(DisorderScore).filter(DisorderScore.isoform_id == isoform_id).first()
        if not score:
            # On-demand calculation if missing
            return self.analyze_isoform(db, isoform_id)
        return score
