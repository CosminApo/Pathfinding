#ifndef _SDL_MANAGER_H
#define _SDL_MANAGER_H

#include "SDL.h" /*for using SDL*/
#include <iostream> /*for inputting and outputting to console*/
#include <stdlib.h>  /*rand, srand*/
#include <vector> /*to store the map*/
#include "SDL_image.h" /*to load images*/


class sdlManager
{
private:
	SDL_Texture* myTex; //the texture of the file
	std::string sprite; //the location of the image 
	SDL_Window* myWindow; //the window to present
	SDL_Renderer* myRenderer; //the rendere used in the window
	bool isRunning; //flag to see if its running
	SDL_Rect myRect; //stores current image
	SDL_Rect myDRect; //stores where to put image

public:
	void init(int xResolution, int yResolution); //function to initialise SDL
	void render(std::vector<std::vector<char>> map, int rows, int columns); //function to render the map
	void createTexture(SDL_Renderer* _renderer, char _currentPos); //function to create the texture
	void clean(); //function to clean all pointers
	bool getisRunning() { return isRunning; }
	sdlManager();
	~sdlManager();
};


#endif // !_GAME_MANAGER_H
