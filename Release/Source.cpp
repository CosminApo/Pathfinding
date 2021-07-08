#include <iostream> /*for inputting and outputting to console*/
#include <random> /*for generating random numbers*/
#include <time.h> /*used to seed the random engine*/
#include <chrono> /*used to seed the random engine*/
#include <cmath> /*used for math functions*/
#include <fstream> /*used to read from file*/
#include <vector> /*used to store data of dynamic size*/
#include <Windows.h> /*used to change the look of the console output*/
#include <algorithm> /*used for vector related functions*/
#include <string> /*used to store chromosomes and other*/
#include <SDL.h> /*used to display the result graphically*/
#include "sdlManager.h" /*the manager for SDL*/

int rows; //the amount of rows in the maze
int columns;  //the amount of columns in the maze
const int xResolution = 1280; //specifies the resolution of the SDL window
const int yResolution = 720;

struct chromosome //structure used to store chromosomes
{
	std::string genes; // string used to store all the genes of the chromosomes (instructions set)
	int xPos; //store the X position of the chromosome in the maze 
	int yPos; //store the Y position of the chromosome in the maze 
	float fitness;  //store the fitness of the chromosome
};

struct node //structure used to store nodes (nodes = grid slots on the maze)
{
	int xPos; //store the X position of the node in the maze
	int yPos; //store the Y position of the node in the maze
	int parent; //store a reference to who is the parent of the node 
	int cost; //store the total cost of the node
	int hcost; //store the H cost of the node

	inline bool operator==(node a) { //overriding the operator "==" for nodes
		if (a.xPos == xPos && a.yPos == yPos) //if the X and Y position of the two nodes being compared matches 
			return true; // then the nodes are equal
		else
			return false; // else they are not equal
	}
};

node findStartingNode(std::vector<std::vector<char>> map) //function used to find the starting node in the maze
{
	node startingNode; //the starting node that will be returned
	for (int i = 0; i < rows; i++) //loop through the maze
	{
		for (int j = 0; j < columns; j++)
		{
			if (map[i][j] == 'S') //find starting node based on char
			{
				startingNode.yPos = i; //set the y position of the node
				startingNode.xPos = j; //set the x position of the node
				startingNode.parent = -1; //original starting point has no parent
				startingNode.cost = 0; //original starting point has no cost
			}
		}
	}
	if (startingNode.parent != -1) //if no node has been found then an error is output
	{
		std::cout << "\n ERROR: No starting node has been found! \n";
	}
	return startingNode;

}

node findEndingNode(std::vector<std::vector<char>> map) //function used to find the ending node in the maze
{
	node endingNode; //the ending node that will be returned
	for (int i = 0; i < rows; i++) //loop through the maze
	{
		for (int j = 0; j < columns; j++)
		{
			if (map[i][j] == 'E') //find starting node based on char
			{
				endingNode.yPos = i; //set the y position of the node
				endingNode.xPos = j; //set the x position of the node
				endingNode.parent = -1; //original ending point has no parent
				endingNode.cost = 0; //original ending point has no cost
			}
		}
	}
	if (endingNode.parent != -1) //if no node has been found then an error is output
	{
		std::cout << "\n ERROR: No ending node has been found! \n";
	}
	return endingNode;

}

node findLowestCost(std::vector<node> &openList) //function used to find the lowest cost node in the open list
{
	int lowestCost = 999; //set to impossible value at start
	std::vector<node> lowestCostNodes; //structure used to store all the lowest cost nodes
	node lowestCostNode; //store the lowest cost node
	for (std::vector<node>::iterator it = openList.begin(); it != openList.end(); ++it) // loop through all nodes in open list
	{
		if (it->cost <= lowestCost) // if the cost of the code examined is lower than the current lowest cost
		{
			lowestCostNodes.push_back(*it); //update the current lowest cost node
			lowestCost = it->cost; //update current lowest cost
		}
	}

	//if there are multiple lowest cost nodes you choose the one with the lowest H cost (closest to ending point)
	int lowestHCost = 999; //set to impossible value at start
	for (std::vector<node>::iterator it = lowestCostNodes.begin(); it != lowestCostNodes.end(); ++it) // loop through all nodes in open list
	{
		if (it->hcost <= lowestHCost) // if the cost of the code examined is lower than the current lowest cost
		{
			lowestCostNode = *it; //update the current lowest cost node
			lowestHCost = it->hcost; //update current lowest cost
		}
	}

	if (lowestCost != 999) //make sure we found a node in the open list
	{
		return lowestCostNode;
	}
	else //if no node has been found output an error message
	{
		std::cout << "ERROR: No nodes weere found in the open list";
	}

}

bool isNodeValid(node testingNode, std::vector<std::vector<char>> map, std::vector<node> closedList, std::vector<node> openList) //function used to find if a node is valid
{
	bool isValid = false; //flag to check if the node is valid

	if (testingNode.yPos >= 0 && testingNode.xPos >= 0 && testingNode.yPos < rows && testingNode.xPos < columns) // make sure node exists within the maze
	{
		if (map[testingNode.yPos][testingNode.xPos] != 'X') //if the node is not an obstacle
		{
			if (std::find(closedList.begin(), closedList.end(), testingNode) == closedList.end()) // check if node is not in the closed list
			{
				if (std::find(openList.begin(), openList.end(), testingNode) == openList.end()) // check if node is not in the open list
				{
					isValid = true; //set the flag to true
				}
			}
		}
	}
	return isValid;
}

std::vector<node> findAdjacentNodes(std::vector<std::vector<char>> map, node currentNode, std::vector<node> closedList, std::vector<node> openList) //function used to find adjacent nodes
{
	std::vector<node> adjacentNodes; //structure to store the adjacent nodes
	node testingNode; //store the current node being tested

// TO GET ALL THE ADJACENT NODES THE X AND Y VALUES SHOULD CHANGE AS FOLLOWS
//		- 1, -1	   - 1, +0     - 1, +1
//
//		+ 0, -1		Y  	X	   - 0, +1
//
//		+ 1, -1    + 1, +0     + 1, +1

	testingNode = currentNode; //set the testing node to the node we're trying to get adjacent nodes for
	testingNode.yPos -= 1; //set the new Y position
	testingNode.xPos -= 1; //set the new X position
	if (isNodeValid(testingNode, map, closedList, openList)) { adjacentNodes.push_back(testingNode); } //add the node to the list if it is valid

	testingNode = currentNode; //reset the testing node to the node we're trying to get adjacent nodes for
	testingNode.yPos -= 1; //set the new Y position
	if (isNodeValid(testingNode, map, closedList, openList)) { adjacentNodes.push_back(testingNode); } //add the node to the list if it is valid

	testingNode = currentNode;
	testingNode.yPos -= 1; //set the new Y position
	testingNode.xPos += 1; //set the new X position
	if (isNodeValid(testingNode, map, closedList, openList)) { adjacentNodes.push_back(testingNode); } //add the node to the list if it is valid

	testingNode = currentNode;
	testingNode.xPos -= 1; //set the new X position
	if (isNodeValid(testingNode, map, closedList, openList)) { adjacentNodes.push_back(testingNode); } //add the node to the list if it is valid

	testingNode = currentNode;
	testingNode.xPos += 1; //set the new X position
	if (isNodeValid(testingNode, map, closedList, openList)) { adjacentNodes.push_back(testingNode); } //add the node to the list if it is valid

	testingNode = currentNode;
	testingNode.yPos += 1; //set the new Y position
	testingNode.xPos -= 1; //set the new X position
	if (isNodeValid(testingNode, map, closedList, openList)) { adjacentNodes.push_back(testingNode); } //add the node to the list if it is valid

	testingNode = currentNode;
	testingNode.yPos += 1; //set the new Y position
	if (isNodeValid(testingNode, map, closedList, openList)) { adjacentNodes.push_back(testingNode); } //add the node to the list if it is valid

	testingNode = currentNode;
	testingNode.yPos += 1; //set the new Y position
	testingNode.xPos += 1; //set the new X position
	if (isNodeValid(testingNode, map, closedList, openList)) { adjacentNodes.push_back(testingNode); } //add the node to the list if it is valid



	return adjacentNodes;

}

int calculateCost(std::vector<node> &adjacentNodes, node startingNode, node endingNode) //function used to calculate the cost of a node
{
	for (std::vector<node>::iterator it = adjacentNodes.begin(); it != adjacentNodes.end(); ++it) // loop through all nodes in the adjacent nodes list
	{
		/*C = S + H;
		Cost = Distance from starting point + Distance from ending point;
		Cost = (Current point – starting point) + (End point – current point).
		S cost = realistic distance while H cost = direct distance (diagonal)
		*/

		int s = abs(it->xPos-startingNode.xPos)+ abs(it->yPos-startingNode.yPos); //find the S cost
		int h = sqrt(pow(abs(endingNode.xPos - it->xPos),2)) + sqrt(pow(abs(endingNode.yPos - it->yPos),2)); //find the H cost
		it->cost = s + h;  // add the two costs together and store them in the chromosome
		it->hcost = h; // store the H cost in the chromosome
	}
	return 0;
}


int aStarSolver(std::vector<std::vector<char>> &map) //function used to solve the maze using A* algorithm
{
	node startingNode; //store the starting node
	node endingNode;  //store the ending node
	std::vector<node> openList;  //structure used to store all the nodes in the open list
	std::vector<node> closedList; //structure used to store all the nodes in the closed list
	std::vector<node> adjacentNodes; //structure used to store all the adjacent nodes
	node currentNode; //store the current node
	int i = 0; //used to keep track of the parent ID
	startingNode = findStartingNode(map); //find stating node
	endingNode = findEndingNode(map); //find ending node
	openList.push_back(startingNode); // add starting node to open list
	while (!openList.empty()) //while the open list is not empty
	{
		adjacentNodes.clear(); //reset the adjacent nodes list
		currentNode = findLowestCost(openList); //find the lowest cost node from the list


		if (currentNode == endingNode) //check if nodes are the same one
		{
			std::cout << "Path complete" << std::endl; //if they are, a path has been found and appropriate message is output
			closedList.push_back(currentNode); //store the current node in the closed list
			break;
		}

		else
		{
			closedList.push_back(currentNode); //add the current node to the closed list
			openList.erase(std::remove(openList.begin(), openList.end(), currentNode), openList.end()); // remove the current node from the open list

			adjacentNodes = findAdjacentNodes(map, currentNode,closedList,openList); //find all the adjacent nodes to the current node

			//SET PARENT NODE
			for (std::vector<node>::iterator it = adjacentNodes.begin(); it != adjacentNodes.end(); ++it) // loop through all nodes in the adjacent nodes list
			{
				it->parent = i; //set the parent's ID of the node
					
			}
			i++; //increase the parent ID
			calculateCost(adjacentNodes, startingNode, endingNode); // calculate the costs of all the adjacent nodes

			// MOVE TO OPEN LIST
			for (std::vector<node>::iterator it = adjacentNodes.begin(); it != adjacentNodes.end(); ++it) // loop through all nodes in the adjacent nodes list
			{
				openList.push_back(*it); //push back the adjacent nodes into the open list
			}

		}
		

	}
	while (currentNode.parent >= -1) //while the current node still has a parent
	{
		if (map[currentNode.yPos][currentNode.xPos] != 'S' && map[currentNode.yPos][currentNode.xPos] != 'E') //if you're not overrading the start or ending operator
		{
			map[currentNode.yPos][currentNode.xPos] = 'V'; //replace the value at map[x][y] with a V to indicate the path found
		}
		
		std::cout << currentNode.xPos << currentNode.yPos << std::endl; //output the current node X and Y positions to show the complete path
		if (currentNode.parent == -1) { break; } //ensure that the code doesnt go forward if the node has no parent 
		currentNode = closedList[currentNode.parent]; //set the current node to its parent
	}

	return 0;
}

std::vector<chromosome> generateChromosomes() //function used to generate chromosomes
{
	std::vector<chromosome> chromosomes; //structure used to store the generated chromosomes
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); //create a seed for the random engine based on current time
	std::default_random_engine generator(seed); //create and seed a random engine
	std::uniform_int_distribution<int> d{ 0, 1 }; //create an uniform int distribution between 0 and 1 that will be used to generate genes
	std::string genes = ""; //string to store the genes generated
	chromosome chrom; //store the current chromosome generated

	for (int i = 0; i < 100; i++) //generating 100 chromosomes  
	{
		genes = ""; //reset chromosome
		for (int j = 0; j < 50; j++) //each chromosome has length 50 (25 possible moves)
		{
			int tempRandNo = d(generator);  //generate either 0 or 1
			genes = genes + std::to_string(tempRandNo); //add bit to gene
		}
		chrom.genes = genes; //add gene to chromosome
		chrom.fitness = 1; //set fitness to max by default
		chromosomes.push_back(chrom); //add the chromosome to the list of chromosomes created
	}
	return chromosomes;
}

bool runChromosome(node startingNode, node endingNode, chromosome &chrom, std::vector<std::vector<char>> &map) //function used to execute the instructions stored in a chromosome
{
	/*
	Instructions encoded per gene:
	00 Up
	01 Right
	10 Down
	11 Left
	*/

	bool solved = false;//flag used to check if the chromosome has solved the solution (found a valid poath)

	chrom.xPos = startingNode.xPos; //reset the X position of the chromosome to the starting position
	chrom.yPos = startingNode.yPos; //reset the Y position of the chromosome to the starting position
	int cposition = 0; //position within chromosome

	while (cposition < chrom.genes.length() - 1) //run through chromosme
	{
		if (chrom.genes[cposition] == '0' && chrom.genes[cposition + 1] == '0') //if genes are 00
		{
			//move up -> decrease the Y position
			int tempX = chrom.xPos; 
			int tempY = chrom.yPos - 1; 

			if (!(tempX < 0 || tempY < 0 || tempX >= columns || tempY >= rows)) //check if move doesn't result in the chromosome being outside the maze
			{
				if (!(map[tempY][tempX] == 'X')) //check if the chromosome doesnt colide with a wall
				{
					if (endingNode.xPos == tempX && endingNode.yPos == tempY) //if you reach the end
					{
						chrom.genes.resize(cposition + 2); //trim the chromosome to remove usless genes
						solved = true; //set the solved flag to true
						break; //end loop
					}
					else //if you dont reach the end
					{
						map[tempY][tempX] = 'V'; //set the value of the map at the current coordinates to track the path
						chrom.xPos = tempX; //set the new X position of the chromosome
						chrom.yPos = tempY; //set the new Y position of the chromosome
					}
				}
			}
		}

		else if (chrom.genes[cposition] == '0' && chrom.genes[cposition + 1] == '1') //if genes are 01
		{
			//move right -> increase the X position
			int tempX = chrom.xPos + 1; 
			int tempY = chrom.yPos;

			if (!(tempX < 0 || tempY < 0 || tempX >= columns || tempY >= rows)) //check if move doesn't result in the chromosome being outside the maze
			{
				if (!(map[tempY][tempX] == 'X')) //check if the chromosome doesnt colide with a wall 
				{
					if (endingNode.xPos == tempX && endingNode.yPos == tempY) //if you reach the end
					{
						chrom.genes.resize(cposition + 2); //trim the chromosome to remove usless genes
						solved = true;//set the solved flag to true
						break; //end loop
					}
					else
					{
						map[tempY][tempX] = 'V'; //set the value of the map at the current coordinates to track the path
						chrom.xPos = tempX; //set the new X position of the chromosomet
						chrom.yPos = tempY; //set the new Y position of the chromosomet
					}
				}
			}
		}

		else if (chrom.genes[cposition] == '1' && chrom.genes[cposition + 1] == '0') //if genes are 10
		{
			//move down -> increase the Y position
			int tempX = chrom.xPos;
			int tempY = chrom.yPos + 1; 

			if (!(tempX < 0 || tempY < 0 || tempX >= columns || tempY >= rows)) //check if move doesn't result in the chromosome being outside the maze
			{
				if (!(map[tempY][tempX] == 'X')) //check if the chromosome doesnt colide with a wall
				{
					if (endingNode.xPos == tempX && endingNode.yPos == tempY) //if you reach the end
					{
						chrom.genes.resize(cposition + 2); //trim the chromosome to remove usless genes
						solved = true; //set the solved flag to true
						break; //end loop
					}
					else
					{
						map[tempY][tempX] = 'V';//set the value of the map at the current coordinates to track the path
						chrom.xPos = tempX;//set the new X position of the chromosome
						chrom.yPos = tempY;//set the new Y position of the chromosome
					}
				}
			}
		}

		else if (chrom.genes[cposition] == '1' && chrom.genes[cposition + 1] == '1') //if genes are 11
		{
			//move left -> decrease the Y position
			int tempX = chrom.xPos - 1; 
			int tempY = chrom.yPos;

			if (!(tempX < 0 || tempY < 0 || tempX >= columns || tempY >= rows)) //check if move doesn't result in the chromosome being outside the maze
			{
				if (!(map[tempY][tempX] == 'X')) //check if the chromosome doesnt colide with a wall
				{
					if (endingNode.xPos == tempX && endingNode.yPos == tempY) //if you reach the end
					{
						chrom.genes.resize(cposition + 2); //trim the chromosome to remove usless genes
						solved = true; //set the solved flag to true
						break; //end loop
					}
					else
					{
						map[tempY][tempX] = 'V'; //set the value of the map at the current coordinates to track the path
						chrom.xPos = tempX;//set the new X position of the chromosome
						chrom.yPos = tempY;//set the new Y position of the chromosome
					}
				}
			}
		}
		cposition += 2; //increase the position within the chromosome so the next instruction is checked
	}
	return solved;
}

int calculateFitness(chromosome &chrom, node endingNode) //function used to calculate the fitness of a chromosome
{
	/*
	F = 1 / (1+dx+dy)
	Fitness = 1 / (1+ the difference in x + the difference in y) 
	Note that we need the absolute distance -> difference between 5 (current pos) and 3 (end pos) = 3-5 = -2 = 2 nodes away
	*/
	int dx = abs(endingNode.xPos - chrom.xPos); //distance from the ending node to current chromosome
	int dy = abs(endingNode.yPos - chrom.yPos); //distance from the ending node to current chromosome
	chrom.fitness = 1.0f / (dx + dy + 1.0f); //set the fitness of the chromosome to the new value 
	return 0;
}
std::vector<chromosome> breedChromosomes(std::vector<chromosome> &chromosomes) //function used to breed chromomse
{
	//GET TOTAL FITNESS
	float totalFitness = 0; //store the total fitness (used to for the roulette wheel method)
	const float crossoverRate = 0.7f; //set the crossover rate (determines if the two chromosomes will swap genes or not)
	const float mutationRate = 0.1f; //set the mutation rate (determines if the two chromosomes's genes will mutate or not)
	std::vector<chromosome> newChroms;
	for (std::vector<chromosome>::iterator it = chromosomes.begin(); it != chromosomes.end(); ++it) // loop through all chromosomes in the list
	{
		totalFitness += (it->fitness); //add up total fitness 
	}

	//SETUP RANDOM
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); //create a seed for the random engine based on current time 
	std::default_random_engine generator(seed); //create and seed a random engine
	std::uniform_real_distribution<float> d{ 0, 1 }; //create an uniform float distribution between 0 and 1 that will be used for the roulette wheel parent selection method
	std::uniform_real_distribution<float> crossover{ 0, 1 }; //create an uniform float distribution between 0 and 1 that will be used to crossover chromosomes
	std::uniform_real_distribution<float> mutate{ 0, 1 }; //create an uniform float distribution between 0 and 1 that will be used to to mutate chromosomes


	for (int i = 0; i < chromosomes.size() / 2; i++) //loop once for every two chromosomes (since 2 chromosomes create 2 children)
	{
		//CHOOSE PARENTS
		chromosome parentOne; //store the parent chromosome 1
		chromosome parentTwo; //store the parent chromosome 2
		chromosome childOne;  //store the child chromosome 1
		chromosome childTwo; //store the child chromosome 2

		float currentPercent = 0; //current percent used for the roulette wheel
		float randomChromosome = d(generator); //generate a random value for the roulette wheel (this determines the chromosome chosen)


		for (std::vector<chromosome>::iterator it = chromosomes.begin(); it != chromosomes.end(); ++it) // loop through all chromosomes in the list
		{
			if (randomChromosome >= currentPercent && randomChromosome <= (it->fitness / totalFitness) + currentPercent) //if the chromosome is chosen
			{
				parentOne = *it; //set the parent 1 chromosome to the current chromosome analysed
				break;
			}
			currentPercent += it->fitness / totalFitness; //add the current fitness to the current percent 
		}

		while (parentTwo.genes == "") //while we dont have a second parent
		{
			currentPercent = 0; //reset value
			randomChromosome = d(generator); //generate a random value for the roulette wheel
			for (std::vector<chromosome>::iterator it = chromosomes.begin(); it != chromosomes.end(); ++it) // loop through all chromosomes in the list
			{
				if (randomChromosome >= currentPercent && randomChromosome <= it->fitness / totalFitness + currentPercent) //if the chromosome is chosen
				{
					if (parentOne.genes != it->genes) //if the chromosome is not equal to the first parent (no asexual reproduction)
					{
						parentTwo = *it; //set the parent 2 chromosome to the current chromosome analysed
						break;
					}
				}
				currentPercent += it->fitness / totalFitness; //add the current fitness to the current percent 
			}
		}

		//CROSSOVER
		std::uniform_int_distribution<int> randomPointCrossove{ 0, 50 }; //create an uniform int distribution between 0 and 1 used to determine the crrossover point
		float randCrossover = crossover(generator); //generate a random number between 0 and 1
		if (randCrossover > 0 && randCrossover < crossoverRate) // determine if a crossover is happening 
		{
			childOne.genes = parentTwo.genes; //set the child 1 genes to the parent 2's genes (swap them around)
			childTwo.genes = parentOne.genes; //set the child 2 genes to the parent 1's genes (swap them around)

			//random point crossover for child one
			for (int i = randomPointCrossove(generator); i < parentOne.genes.length(); i++) //loop between the random point and the length of the chromosome
			{
				childOne.genes[i] = parentOne.genes[i]; //crossover genes
			}
			//random point crossover for child two
			for (int i = randomPointCrossove(generator); i < parentTwo.genes.length(); i++) //loop between the random point and the length of the chromosome
			{

				childTwo.genes[i] = parentTwo.genes[i]; //crossover genes
			}
		}
		else //if no crossover happened
		{
			childOne.genes = parentTwo.genes;
			childTwo.genes = parentOne.genes;
		}

		//MUTATE

		//child one
		for (int i = 0; i < childOne.genes.length(); i++) //loop through the chromosome
		{
			float randMutate = crossover(generator); //generate random number between 0 and 1
			if (randMutate > 0 && randMutate < mutationRate) //if a mutation happens
			{
				//swap genes around
				if (childOne.genes[i] == '1') { childOne.genes[i] = '0'; }
				if (childOne.genes[i] == '0') { childOne.genes[i] = '1'; }

			}
		}
		//child two
		for (int i = 0; i < childTwo.genes.length(); i++)
		{
			float randMutate = crossover(generator); //generate random number between 0 and 1
			if (randMutate > 0 && randMutate < mutationRate) //if a mutation happens
			{
				//swap genes around
				if (childTwo.genes[i] == '1') { childTwo.genes[i] = '0'; }
				if (childTwo.genes[i] == '0') { childTwo.genes[i] = '1'; }

			}
		}

		//move new chromosomes in the list
		newChroms.push_back(childOne);
		newChroms.push_back(childTwo);
	}


	return newChroms;
}

int GAsolver(std::vector<std::vector<char>> &map) //function used to solve the maze using Genetic Algorithm
{
	std::vector<chromosome>  chromosomes; //structure to store the chromosomes
	chromosomes = generateChromosomes(); //generate the first batch of chromosomes
	node endingNode = findEndingNode(map); //find the ending node
	node startingNode = findStartingNode(map); //find the starting node
	bool solved = false; //flag used to check if a solution has been reached
	chromosome finalChromosome; //store the final chromosome (the one that contains a solution)
	std::vector<std::vector<char>> secondMap = map; //structure used to store a copy of the map (to ensure that only the path of the final chromosome is highlighted)
	int i = 0; //variable used to indicate the current chromosome generation

	while (!solved)
	{
		secondMap = map; //reset the map copy
		std::cout << "Gen " << i << " chromosomes: " << chromosomes.size() << std::endl; //output the current generation of chromosomes
		i++; //increase the generation number
		for (std::vector<chromosome>::iterator it = chromosomes.begin(); it != chromosomes.end(); ++it) // loop through all chromosomes in the list
		{
			solved = runChromosome(startingNode, endingNode, *it, secondMap); //run each chromosome
			if (solved) //if a solution has been reached
			{
				finalChromosome = *it; //set the final chromosome to the current chromosome analysed
				break;
			}
			calculateFitness(*it, endingNode); //get the fitness for each chromosome
		}
		if (solved) //make sure the code doesnt advance if a solution has been found (failsafe)
		{
			break;
		}
		float fitnessTemp = 0; //store the a temporary fitness to check which chromosome has the highest fitness
		for (std::vector<chromosome>::iterator it = chromosomes.begin(); it != chromosomes.end(); ++it) // loop through all chromosomes in the list
		{
			if (it->fitness > fitnessTemp) //if the current chromosome's fitness is bigger than the current biggest fitness
			{
				finalChromosome = *it; //set the best chromosome to the current chromosome
				fitnessTemp = it->fitness; //set the new fitness value
			}
		}
		std::cout << " Best Chromosome: " << finalChromosome.genes << " Fitness: " << fitnessTemp << std::endl; // output the best chromosome
		std::vector<chromosome> temp = breedChromosomes(chromosomes); //generate 20 new chromosomes
		chromosomes.clear(); //clear the current list of chromosomes
		chromosomes = temp; //set the new chromosomes to the list 

	}
	std::cout << finalChromosome.genes; //output the solution
	runChromosome(startingNode, endingNode, finalChromosome, map); //run the chromosome to highlight the path

	return 0;
}


std::vector<std::vector<char>> mapCreator(char all[], int size) //function used to create the map from the file
{
	//GET ROWS AND COLUMNS
	std::string tempRows; //store a temporary value for rows
	std::string tempColumns; //store a temporary value for chroms
	

	int i = 0, j = 0;  //create two variables used to sanitize the input and find rows and columns

	while (!(all[i + 1] == '\0' && all[i + 2] == '\0')) //read in the current columns 
	{
		if (all[i] > 47 && all[i] < 58) //ensure its a numeric character
		{
			tempColumns += all[i];
			j++;
		}
		i++;
	}
	tempColumns += all[i]; 

	i++;
	while (!(all[i + 1] == '\0' && all[i + 2] == '\0')) //read in the current rows
	{
		if (all[i] > 47 && all[i] < 58) //ensure its a numeric character
		{
			tempRows += all[i];
		}
		i++;
	}
	tempRows += all[i]; 

	rows = std::stoi(tempRows); //set the global rows value to the temp one
	columns = std::stoi(tempColumns); //set the global columns value to the temp one

	for (int i = 0 , j = 0; i <= size - tempColumns.size() - tempRows.size(); i++) //loop through the input to remove all symbols
	{
		if (all[i] > 47 && all[i] < 58) //if the char is not a number
		{
			all[j] = all[i]; //remove it
			j++;
		}
	}
	//

	int index = (tempColumns.size() + tempRows.size()); // start at 2 to avoide size defining chars
	std::vector<std::vector<char>> map(rows, std::vector<char>(columns, 0)); //create 2d vector for map
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			if (all[index] == '0')
			{
				map[i][j] = '0'; // walkable 
			}
			else if (all[index] == '1')
			{
				map[i][j] = 'X'; //wall
			}
			else if (all[index] == '2')
			{
				map[i][j] = 'S'; //start
			}
			else if (all[index] == '3')
			{
				map[i][j] = 'E'; //end
			}

			index++;
		}
	}
	return map;
}
int mapDisplayer(std::vector<std::vector<char>> map) //function used to display the map
{ 
	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get a handle for the console (used to change colour)
	for (int i = 0; i < rows; i++) //loop through the map
	{	
		for (int j = 0; j < columns; j++)
		{
			if (map[i][j] == 'V') //if it is part of the Path
			{
				SetConsoleTextAttribute(hConsole, 10); //change the colour of the console text to green for better visibility
			}
			std::cout << map[i][j]; // output from map
			SetConsoleTextAttribute(hConsole, 7); //reset teh colour of the console
		}
		std::cout << std::endl; //separate rows
	}
	return 0;
}

int main(int argc, char *argv[]) //main function 
{
	sdlManager manager;
	manager.init(xResolution, yResolution); //initialise the game
	std::ifstream inFile; //store the current file here
	char x; //store the currently read char from the file
	char all[500]; //store all the chars from the file
	std::string filePath; //store the path of the file
	char algChoice; //store which algorithm to use
	int size = 0; //store the size of the input

	while (size <= 0) //run untill it finds a file
	{
		std::cout << "Please insert the name of the file containing the maze (including .txt): " << std::endl;
		std::cin >> filePath; //store user input

		inFile.open(filePath); //open the file
		while (inFile >> x) //read in characters while possible
		{
			all[size] = x; //store the character 
			size++;

		}
		if (size == 0)
		{
			std::cout << "\nERROR: file not found\n";
		}
	}

	std::cout << "Please select Genetic Algorithm (1) or A* pathfinding (2) " << std::endl; 
	std::cin >> algChoice; //store user input
	
	const int constSize = size; 
	std::vector<std::vector<char>> map; //store the map 
	map = mapCreator(all,size); //generate the map
	
	mapDisplayer(map); //display the map in the console
	
	std::cout << std::endl;
	std::cout << std::endl;

	if (algChoice == '1') //choose which algorithm to use
	{
		GAsolver(map);
	}
	else
	{
		aStarSolver(map);
	}
	
	std::cout << std::endl;
	std::cout << std::endl;

	mapDisplayer(map); //display the solved map in the console

	while (manager.getisRunning()) //keep the SDL window open
	{
		manager.render(map, rows, columns); //render the map
	}

	manager.clean(); //clean SDL elements
 	inFile.close(); //close file
	return 0;

}