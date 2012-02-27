
#include"PressTestEngine.hpp"
#include<yade/core/Scene.hpp>


void PressTestEngine::action(){
	if (curentDirection != stop) {
		if (curentDirection==forward) { 										 ///<Forward direction of the press
			FOREACH(Body::id_t id, ids){
				if (!(scene->bodies->exists(id))) continue;
				currentVerticalForce = scene->forces.getForce(id)[2]; 	///<Define current vertical force
				minimalForce = maxVerticalForce*0.1;									///<Define minimal edge of the force (10% from Maximal)
				minimalPredictedForce = predictedForce*0.1;						///<Define minimal edge of the Predicted force (10% from Predicted)
				if (currentVerticalForce > maxVerticalForce) {				///<Force increasing. Press is working normally
					maxVerticalForce = currentVerticalForce;
				} else if ((currentVerticalForce<=(minimalForce))&&(maxVerticalForce>minimalPredictedForce)) {
					/**
					 * Force is decreased lower, than minimal. The body seems "cracked".
					 * Starting the countdown
					 */ 
					currentIterationAfterDestruction++;
					if (currentIterationAfterDestruction>=numberIterationAfterDestruction) {
						/**
						 * The body definitly cracked. Change the direction of press and increase the velosity in 5 times.
						 */ 
						curentDirection=backward;
						TranslationEngine::velocity *= -pressVelocityForw2Back;
						currentIterationAfterDestruction = scene->iter/pressVelocityForw2Back*riseUpPressHigher;
					}
				}  else if (((currentIterationAfterDestruction!=0)&&(maxVerticalForce !=0))) {
					/**
					 * We have found, that it was not "Finish destruction"
					 * Abnulling currentIterationAfterDestruction
					 */
					currentIterationAfterDestruction=0;
				}
			}
			TranslationEngine::action();
		} else if (curentDirection==backward) {							 ///<The press returns back to the normal position
			if (currentIterationAfterDestruction > 0) {
				currentIterationAfterDestruction--;
				TranslationEngine::action();
			} else {
				curentDirection=stop;														///<If the press is in normal position -> STOP
				Omega::instance().pause();					///<Stop simulation loop
			}
		}
	}
}

CREATE_LOGGER(PressTestEngine);
YADE_PLUGIN((PressTestEngine));

