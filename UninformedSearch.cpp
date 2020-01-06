// Code for Project 1 of CSCI404A
// Utkarsh Srivastava (10832610)

// Defining Header file

#include <stdio.h>
#include <vector>
#include <fstream>
#include <iterator>
#include <queue>
#include <limits>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <string>

using namespace std;

// Initialization

int main(int argc, char** argv){

    string start_state = argv[2];
    string goal_state = argv[3];
    vector<string> city_list;
    vector< pair<int,int> > city_dist[200];
    int split_counter = 0;
    int item_counter = 0;
    
    string line;

// taking two temperory variable temp1 and temp2 having distance as dist1 and dist2

    string split;
    string temp1, temp2;
    vector<string>::iterator temp_city1, temp_city2;
    int dist, dist1, dist2;

    ifstream map(argv[1]);
    if (map.is_open()){
        while ( getline (map,line) ){
                stringstream s2ss(line);
                
                if (!line.compare("END OF INPUT") ){
                    break;
                }

                while( getline (s2ss,split,' ') ){
                    
                    if(split_counter%3 == 0){
                        temp1 = split;
                        
                    }
                    else if(split_counter%3 == 1){
                        temp2 = split;
                        
                    }

                    if(split_counter%3 < 2){
                        if(find(city_list.begin(),city_list.end(),split) == city_list.end() ){
                            city_list.push_back(split);     
                        }
                    }
                    else{

                        // got temporary value for the city and the distance between them

                        temp_city1 = find(city_list.begin(),city_list.end(),temp1);
                        temp_city2 = find(city_list.begin(),city_list.end(),temp2);

                        stringstream s2num(split);
                        s2num >> dist;

                        dist1 = std::distance(city_list.begin(),temp_city1);
                        dist2 = std::distance(city_list.begin(),temp_city2);
                        
                        city_dist[dist1].push_back(make_pair(dist2,dist));
                        city_dist[dist2].push_back(make_pair(dist1,dist));                  
                    }
                    split_counter++;
                }
                split_counter = 0;
                item_counter++;
        }
        map.close();
    }

// check for start and goal state i.e cities and if no name is found then return "City name not found in input file, check city name"


  int start_state_it,goal_state_it;
  if(find(city_list.begin(), city_list.end(), start_state) != city_list.end() && find(city_list.begin(), city_list.end(), goal_state) != city_list.end()){
        start_state_it = distance(city_list.begin(),find(city_list.begin(), city_list.end(), start_state));
        goal_state_it = distance(city_list.begin(),find(city_list.begin(), city_list.end(), goal_state));
  }
  else{
        cout << "City name not found in input file, check city name" <<'\n';
        return 0;
  }
  vector <pair <int, int> > tempvec[city_list.size()];

    for(int j = 0; j<city_list.size(); j++){
        tempvec[j] = city_dist[j];
  }

// output  the start and goal city

    //cout << "Start State: "<< start_state_it << " Goal State: " << goal_state_it << '\n';
    
    priority_queue<pair<int,int>, vector<pair<int, int> >, greater<pair<int,int> > >  Frontier; 

    Frontier.push(make_pair(0, start_state_it));

    int node, weight, int_state;

    vector<int>path(city_list.size(), std::numeric_limits<int>::max());
    path[start_state_it] = 0;  

    vector<pair <int, int> > answervec1;
    vector<int> answervec2, answervec2w;

// algorithm for finding the distance

    while(!Frontier.empty()){
        int current_state = Frontier.top().second;
        Frontier.pop();
        int node, weight;
        
        for(vector < pair<int, int> >::const_iterator it = tempvec[current_state].begin(); it!=tempvec[current_state].end();it++){
            node = it->first;
            weight = it->second;

            if(path[node] > path[current_state] + weight){
                
                path[node] = path[current_state]+weight;

                answervec1.push_back(make_pair(current_state, node));
                

                if(node == goal_state_it){
                    int_state = goal_state_it;
                
                    int last,next;

                    while(int_state != start_state_it){

                        for(vector < pair<int, int> >::const_iterator iter = answervec1.begin(); iter!=answervec1.end();iter++){
                            last = iter->first;
                            next = iter->second;

                            if(int_state == next){
                                answervec2.push_back(int_state);
                                
                                for(vector < pair<int, int> >::const_iterator d = tempvec[last].begin(); d!=tempvec[last].end();d++){
                                    int oldn = d->first;
                                    int oldw = d->second;
                                    if (oldn == int_state){
                                        
                                        answervec2w.push_back(oldw);    
                                    }
                                }
                                int_state = last;
                            }
                        }
                    }
                    
                }
                Frontier.push(make_pair(path[node], node));

            
            }

    
        }
    }
for (int k = 0; k < city_list.size(); k++){

}

answervec2.push_back(start_state_it);

// Displaying the route:

if(path[goal_state_it] == std::numeric_limits<int>::max()){
    cout << "Distance: Infinity" <<'\n';
    cout << "Route: "<< '\n';
    cout << "None "<<'\n';
}

else{
cout << "Distance: " << path[goal_state_it] << "km" <<'\n';
cout << "Route: "<< '\n';
}


for(int w = answervec2.size()-1; w > 0;w--){
            
    cout << city_list[answervec2[w]] << " to "<< city_list[answervec2[w-1]] << ", " <<answervec2w[w-1] <<" km" << '\n';
}
    
return 0;
}