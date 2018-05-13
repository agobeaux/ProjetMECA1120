/**
	Fonction qui renvoie soit l'élément contenant la balle située en (x,y)
	soit, le voisin situé entre l'ancien élément la contenant et la position
	de la balle actuelle
*/
int elemContains(double x, double y, femMesh *theMesh, int iElem){
	double jacobian[3];
	double x2, x3, y2, y3;
	int tab2[3] = {0,0,1};
	int tab3[3] = {1,2,2};
	for(int i = 0; i < 3; i++){
		x2 = theMesh->X[3*iElem+tab2[i]];
		y2 = theMesh->Y[3*iElem+tab2[i]];
		x3 = theMesh->X[3*iElem+tab3[i]];
		y3 = theMesh->Y[3*iElem+tab3[i]];
		jacobian[i] = ((x - x3) * (y2 - y3) - (x2 - x3) * (y - y3) >= 0);
	}
	if(jacobian[0] == jacobian[1] && jacobian[1] == jacobian[2]){
		return iElem;
	}
	else{
		if(jacobian[0] == jacobian[1]){ // jacobian[2] d'un autre signe
			return theMesh->neighbours[3*iElem+2];
		}
		elseif(jacobian[1] == jacobian[2]){ // jacobian[0] d'un autre signe
			return theMesh->neighbours[3*iElem+0];
		}
		else{ // jacobian[1] d'un autre signe
			return theMesh->neighbours[3*iElem+1];
		}
	}
}
