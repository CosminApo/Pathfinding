#include "sdlManager.h"



sdlManager::sdlManager() //constructor
{
}


sdlManager::~sdlManager() //deconstructor
{
}
void sdlManager::init(int xResolution, int yResolution) //function to initialise SDL
{
	isRunning = true; //set the flag to true
	myWindow = nullptr;
	myRenderer = nullptr;
	int initOk = SDL_Init(SDL_INIT_EVERYTHING); //initialise SDL
	if (initOk != 0) //check if SDL has been initialised
	{
		std::cout << "Unable to initialise SDL: " << SDL_GetError() << std::endl;
		isRunning = false;
	}
	else
	{
		std::cout << "SDL Initialised" << std::endl;
	}


	myWindow = SDL_CreateWindow( //init Window
		"Main", //name of the window
		SDL_WINDOWPOS_UNDEFINED, //position of the window
		SDL_WINDOWPOS_UNDEFINED,
		xResolution, yResolution, //size
		SDL_WINDOW_RESIZABLE);

	if (myWindow == NULL) //check if window has been initialised
	{
		std::cout << "Could not create window: " << SDL_GetError << std::endl;
		isRunning = false;
	}
	else
	{
		std::cout << "Window Created" << std::endl;
	}

	myRenderer = SDL_CreateRenderer(myWindow, -1, SDL_RENDERER_ACCELERATED); //initialise Renderer
	if (myRenderer == NULL) //check if renderer has been initialised
	{
		std::cout << "Could not create renderer: " << SDL_GetError << std::endl;
		isRunning = false;
	}
	else
	{
		std::cout << "Renderer Created" << std::endl;
	}

	myRect.w = 32; //set the starting value of the rectangles
	myRect.h = 32;
	myDRect.w = 32;
	myDRect.h = 32;


	SDL_UpdateWindowSurface(myWindow);
}

void  sdlManager::createTexture(SDL_Renderer* _renderer, char _currentPos) //function to create the texture
{ 
	//BASED ON THE CHARACTER ANALYSED ASSIGN A DIFFERENT IMAGE
	if (_currentPos == '0') 
	{
		sprite = "Assets/whiteSquare.png";
	}
	else if (_currentPos == 'X')
	{
		sprite = "Assets/blackSquare.png";
	}
	else if (_currentPos == 'S')
	{
		sprite = "Assets/blueSquare.png";
	}
	else if (_currentPos == 'E')
	{
		sprite = "Assets/redSquare.png";
	}
	else if (_currentPos == 'V')
	{
		sprite = "Assets/greenSquare.png";
	}

	SDL_Surface* tmpsurf = IMG_Load(sprite.c_str()); //create a new surface with the path of the image

	if (tmpsurf == nullptr) { //check if the surface has been initialised
		std::cout << "IMG_Load: " << IMG_GetError() << "\n";
	}
	myTex = SDL_CreateTextureFromSurface(_renderer, tmpsurf); //create texture
	SDL_FreeSurface(tmpsurf);

}

void sdlManager::render(std::vector<std::vector<char>> map, int rows, int columns) //function to render the map
{

	myDRect.x = 0; //set initial values of destination rect
	myDRect.y = 0;
	for (int i = 0; i < rows; i++) //loop through the map
	{
		myDRect.x = 0;
		for (int j = 0; j < columns; j++)
		{
			createTexture(myRenderer, map[i][j]); //create a texture for this position
			SDL_RenderCopyEx(myRenderer, myTex, &myRect, &myDRect, 0, 0, SDL_FLIP_NONE); //copy the texture to the renderer
			myDRect.x += 32; //move the rectangle one column right

			SDL_DestroyTexture(myTex); //clean the  texture
		}
		myDRect.y += 32; //move the rectangle one row down
	}
		

	


	SDL_RenderPresent(myRenderer); //display the renderer
	SDL_RenderClear(myRenderer); //clean the renderer
}

void sdlManager::clean() //function to clean all pointers
{
	SDL_DestroyRenderer(myRenderer);
	SDL_DestroyWindow(myWindow);
	IMG_Quit();
}
