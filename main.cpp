#include "main.h"

//-----------------------------------------------------------------------------
void MyGlDraw(void)
{
	//*************************************************************************
	// Chame aqui as funções do mygl.h
	//*************************************************************************
    PipeLine();
    lerArquivo(&pontos, &faces);
    limpaTela();
    DesenhaMacaco();
}

//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    Face *ptrFace;
    PtrPonto *ptrP;
    for(ptrFace = faces ; ptrFace != NULL ; ptrFace = ptrFace->proximo){
       for(ptrP = ptrFace->ptrP ; ptrP != NULL ; ptrP = ptrP->proximo){
            printf("v  %f  %f  %f\n",ptrP->ponto->x, ptrP->ponto->y, ptrP->ponto->z);
        }
		printf("----------------------------------------------\n");
    }

	// Inicializações.
	InitOpenGL(&argc, argv);
	InitCallBacks();
	InitDataStructures();

	// Ajusta a função que chama as funções do mygl.h
	DrawFunc = MyGlDraw;

	// Framebuffer scan loop.
	glutMainLoop();

	return 0;
}

