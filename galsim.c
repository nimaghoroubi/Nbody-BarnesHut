#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "graphics.h"
#define depth 0


///graphics
int graphics;
const float W = 1;
const float H = 1;
const int   redColor    = 0;
const int   greenColor  = 1;
const int   yellowColor = 2;

typedef struct vector{
    double x;
    double y;
}vec;

typedef struct body_attribute{
    vec position;
    vec velocity;
    double mass;
    double brightness;
}body_attribute;///body struct

typedef struct dummy_attribute{
    vec position;
    double mass;
    int id;
}dummy_attribute;

typedef struct Node_structure{
    vec     quadrant_center;
    vec     quadrant_mass_center;
    double  quadrant_size;
    double  quadrant_mass;
    int     quadrant_body_count;

    struct body_attribute** bodies;

    struct Node_structure* nw;
    struct Node_structure* ne;
    struct Node_structure* sw;
    struct Node_structure* se;

}Node;///node structure





void create_tree(body_attribute** body, Node* parent){

    if (parent->quadrant_body_count == 0){

        parent->bodies = body;


    }else if (parent->quadrant_body_count == 1){
        body_attribute* parent_bodies = *(parent->bodies);
//moving the current body to a sub Node
        if (parent_bodies->position.y <= parent->quadrant_center.y){
            if (parent_bodies->position.x <= parent->quadrant_center.x){
//sw
                Node* node0 = (Node*)calloc(1,sizeof(Node));
                node0->quadrant_size = parent->quadrant_size/2;
                node0->quadrant_center.x = parent->quadrant_center.x - parent->quadrant_size/4;
                node0->quadrant_center.y = parent->quadrant_center.y - parent->quadrant_size/4;
                node0->quadrant_body_count = 1;
                node0->quadrant_mass_center.x = parent_bodies->position.x;
                node0->quadrant_mass_center.y = parent_bodies->position.y;
                node0->quadrant_mass          = parent_bodies->mass;
                node0->bodies = parent->bodies;

                parent->sw = node0;
                parent->bodies = NULL;
                if(depth) DrawRectangle(node0->quadrant_center.x - (node0->quadrant_size)/2, node0->quadrant_center.y - (node0->quadrant_size)/2, W, H, node0->quadrant_size, node0->quadrant_size,0, yellowColor);

            }else{
//se
                Node* node1 = (Node*)calloc(1,sizeof(Node));
                node1->quadrant_size = parent->quadrant_size/2;
                node1->quadrant_center.x = parent->quadrant_center.x + parent->quadrant_size/4;
                node1->quadrant_center.y = parent->quadrant_center.y - parent->quadrant_size/4;
                node1->quadrant_body_count = 1;
                node1->quadrant_mass_center.x = parent_bodies->position.x;
                node1->quadrant_mass_center.y = parent_bodies->position.y;
                node1->quadrant_mass          = parent_bodies->mass;
                node1->bodies = parent->bodies;

                parent->se = node1;
                parent->bodies = NULL;
                if(depth) DrawRectangle(node1->quadrant_center.x - (node1->quadrant_size)/2, node1->quadrant_center.y - (node1->quadrant_size)/2, W, H, node1->quadrant_size, node1->quadrant_size,0, yellowColor);


            }
        }else{
            if (parent_bodies->position.x <= parent->quadrant_center.x){
//nw
                Node* node2 = (Node*)calloc(1,sizeof(Node));
                node2->quadrant_size = parent->quadrant_size/2;
                node2->quadrant_center.x = parent->quadrant_center.x - parent->quadrant_size/4;
                node2->quadrant_center.y = parent->quadrant_center.y + parent->quadrant_size/4;
                node2->quadrant_body_count = 1;
                node2->quadrant_mass_center.x = parent_bodies->position.x;
                node2->quadrant_mass_center.y = parent_bodies->position.y;
                node2->quadrant_mass          = parent_bodies->mass;
                node2->bodies = parent->bodies;

                parent->nw = node2;
                parent->bodies = NULL;
                if(depth) DrawRectangle(node2->quadrant_center.x - (node2->quadrant_size)/2, node2->quadrant_center.y - (node2->quadrant_size)/2, W, H, node2->quadrant_size, node2->quadrant_size,0, yellowColor);

            }else{
//ne
                Node* node3 = (Node*)calloc(1,sizeof(Node));
                node3->quadrant_size = parent->quadrant_size/2;
                node3->quadrant_center.x = parent->quadrant_center.x + parent->quadrant_size/4;
                node3->quadrant_center.y = parent->quadrant_center.y + parent->quadrant_size/4;
                node3->quadrant_body_count = 1;
                node3->quadrant_mass_center.x = parent_bodies->position.x;
                node3->quadrant_mass_center.y = parent_bodies->position.y;
                node3->quadrant_mass          = parent_bodies->mass;
                node3->bodies = parent->bodies;

                parent->ne = node3;
                parent->bodies = NULL;
                if(depth) DrawRectangle(node3->quadrant_center.x - (node3->quadrant_size)/2, node3->quadrant_center.y - (node3->quadrant_size)/2, W, H, node3->quadrant_size, node3->quadrant_size,0, yellowColor);

            }
        }

//moving the new body to a sub Node

        if ((*body)->position.y <= parent->quadrant_center.y){
            if ((*body)->position.x <= parent->quadrant_center.x){
//sw
                if (parent->sw != NULL){
                    create_tree(body, parent->sw);


                }else{
                    Node* node0 = (Node*)calloc(1,sizeof(Node));
                    node0->quadrant_size = parent->quadrant_size/2;
                    node0->quadrant_center.x = parent->quadrant_center.x - parent->quadrant_size/4;
                    node0->quadrant_center.y = parent->quadrant_center.y - parent->quadrant_size/4;
                    node0->quadrant_body_count = 1;
                    node0->quadrant_mass_center.x = (*body)->position.x;
                    node0->quadrant_mass_center.y = (*body)->position.y;
                    node0->quadrant_mass          = (*body)->mass;

                    node0->bodies = body;

                    parent->sw = node0;

                    if(depth) DrawRectangle(node0->quadrant_center.x - (node0->quadrant_size)/2, node0->quadrant_center.y - (node0->quadrant_size)/2, W, H, node0->quadrant_size, node0->quadrant_size,0, yellowColor);

                }
            }else{
//se
                if (parent->se != NULL){
                    create_tree(body, parent->se);

                }else{
                    Node* node1 = (Node*)calloc(1,sizeof(Node));
                    node1->quadrant_size = parent->quadrant_size/2;
                    node1->quadrant_center.x = parent->quadrant_center.x + parent->quadrant_size/4;
                    node1->quadrant_center.y = parent->quadrant_center.y - parent->quadrant_size/4;
                    node1->quadrant_body_count = 1;
                    node1->quadrant_mass_center.x = (*body)->position.x;
                    node1->quadrant_mass_center.y = (*body)->position.y;
                    node1->quadrant_mass          = (*body)->mass;
                    node1->bodies = body;

                    parent->se = node1;

                    if(depth) DrawRectangle(node1->quadrant_center.x - (node1->quadrant_size)/2, node1->quadrant_center.y - (node1->quadrant_size)/2, W, H, node1->quadrant_size, node1->quadrant_size,0, yellowColor);

                }
            }
        }else{
            if ((*body)->position.x <= parent->quadrant_center.x){
//nw
                if (parent->nw != NULL){
                    create_tree(body, parent->nw);


                }else{
                    Node* node2 = (Node*)calloc(1,sizeof(Node));
                    node2->quadrant_size = parent->quadrant_size/2;
                    node2->quadrant_center.x = parent->quadrant_center.x - parent->quadrant_size/4;
                    node2->quadrant_center.y = parent->quadrant_center.y + parent->quadrant_size/4;
                    node2->quadrant_body_count = 1;
                    node2->quadrant_mass_center.x = (*body)->position.x;
                    node2->quadrant_mass_center.y = (*body)->position.y;
                    node2->quadrant_mass          = (*body)->mass;
                    node2->bodies = body;


                    parent->nw = node2;

                    if(depth) DrawRectangle(node2->quadrant_center.x - (node2->quadrant_size)/2, node2->quadrant_center.y - (node2->quadrant_size)/2, W, H, node2->quadrant_size, node2->quadrant_size,0, yellowColor);

                }
            }else{
//ne
                if (parent->ne != NULL){
                    create_tree(body, parent->ne);

                }else{
                    Node* node3 = (Node*)calloc(1,sizeof(Node));
                    node3->quadrant_size = parent->quadrant_size/2;
                    node3->quadrant_center.x = parent->quadrant_center.x + parent->quadrant_size/4;
                    node3->quadrant_center.y = parent->quadrant_center.y + parent->quadrant_size/4;
                    node3->quadrant_body_count = 1;
                    node3->quadrant_mass_center.x = (*body)->position.x;
                    node3->quadrant_mass_center.y = (*body)->position.y;
                    node3->quadrant_mass          = (*body)->mass;
                    node3->bodies = body;


                    parent->ne = node3;

                    if(depth) DrawRectangle(node3->quadrant_center.x - (node3->quadrant_size)/2, node3->quadrant_center.y - (node3->quadrant_size)/2, W, H, node3->quadrant_size, node3->quadrant_size,0, yellowColor);

                }
            }
        }
///mass center of quadrant also, check it


    }else if (parent->quadrant_body_count > 1){

        if ((*body)->position.y <= parent->quadrant_center.y){
            if ((*body)->position.x <= parent->quadrant_center.x){
//sw
                if (parent->sw != NULL){
                    create_tree(body, parent->sw);


                }else{
                    Node* node0 = (Node*)calloc(1,sizeof(Node));
                    node0->quadrant_size = parent->quadrant_size/2;
                    node0->quadrant_center.x = parent->quadrant_center.x - parent->quadrant_size/4;
                    node0->quadrant_center.y = parent->quadrant_center.y - parent->quadrant_size/4;
                    node0->quadrant_body_count = 1;
                    node0->quadrant_mass_center.x = (*body)->position.x;
                    node0->quadrant_mass_center.y = (*body)->position.y;
                    node0->quadrant_mass          = (*body)->mass;
                    node0->bodies = body;


                    parent->sw = node0;

                    if(depth) DrawRectangle(node0->quadrant_center.x - (node0->quadrant_size)/2, node0->quadrant_center.y - (node0->quadrant_size)/2, W, H, node0->quadrant_size, node0->quadrant_size,0, yellowColor);

                }
            }else{
//se
                if (parent->se != NULL){
                    create_tree(body, parent->se);

                }else{
                    Node* node1 = (Node*)calloc(1,sizeof(Node));
                    node1->quadrant_size = parent->quadrant_size/2;
                    node1->quadrant_center.x = parent->quadrant_center.x + parent->quadrant_size/4;
                    node1->quadrant_center.y = parent->quadrant_center.y - parent->quadrant_size/4;
                    node1->quadrant_body_count =1;
                    node1->quadrant_mass_center.x = (*body)->position.x;
                    node1->quadrant_mass_center.y = (*body)->position.y;
                    node1->quadrant_mass          = (*body)->mass;
                    node1->bodies = body;


                    parent->se = node1;
                    if(depth) DrawRectangle(node1->quadrant_center.x - (node1->quadrant_size)/2, node1->quadrant_center.y - (node1->quadrant_size)/2, W, H, node1->quadrant_size, node1->quadrant_size,0, yellowColor);

                }
            }
        }else{
            if ((*body)->position.x <= parent->quadrant_center.x){
//nw
                if (parent->nw != NULL){
                    create_tree(body, parent->nw);

                }else{
                    Node* node2 = (Node*)calloc(1, sizeof(Node));
                    node2->quadrant_size = parent->quadrant_size/2;
                    node2->quadrant_center.x = parent->quadrant_center.x - parent->quadrant_size/4;
                    node2->quadrant_center.y = parent->quadrant_center.y + parent->quadrant_size/4;
                    node2->quadrant_body_count =1;
                    node2->quadrant_mass_center.x = (*body)->position.x;
                    node2->quadrant_mass_center.y = (*body)->position.y;
                    node2->quadrant_mass          = (*body)->mass;
                    node2->bodies = body;


                    parent->nw = node2;

                    if(depth) DrawRectangle(node2->quadrant_center.x - (node2->quadrant_size)/2, node2->quadrant_center.y - (node2->quadrant_size)/2, W, H, node2->quadrant_size, node2->quadrant_size,0, yellowColor);

                }
            }else{
//ne
                if (parent->ne != NULL){
                    create_tree(body, parent->ne);

                }else{
                    Node* node3 = (Node*)calloc(1,sizeof(Node));
                    node3->quadrant_size = parent->quadrant_size/2;
                    node3->quadrant_center.x = parent->quadrant_center.x + parent->quadrant_size/4;
                    node3->quadrant_center.y = parent->quadrant_center.y + parent->quadrant_size/4;
                    node3->quadrant_body_count = 1;
                    node3->quadrant_mass_center.x = (*body)->position.x;
                    node3->quadrant_mass_center.y = (*body)->position.y;
                    node3->quadrant_mass          = (*body)->mass;
                    node3->bodies = body;


                    parent->ne = node3;

                    if(depth) DrawRectangle(node3->quadrant_center.x - (node3->quadrant_size)/2, node3->quadrant_center.y - (node3->quadrant_size)/2, W, H, node3->quadrant_size, node3->quadrant_size,0, yellowColor);

                }
            }
        }

//here

    }


    parent->quadrant_body_count++;
    parent->quadrant_mass += (*body)->mass;
    parent->quadrant_mass_center.x = (((*body)->mass*(*body)->position.x)+
                                      ((parent->quadrant_mass-(*body)->mass)*parent->quadrant_mass_center.x))/parent->quadrant_mass;
    parent->quadrant_mass_center.y = (((*body)->mass*(*body)->position.y)+
                                      ((parent->quadrant_mass-(*body)->mass)*parent->quadrant_mass_center.y))/parent->quadrant_mass;


}

void freetree(Node* root){
    //
    if (root->nw != NULL){
        freetree(root->nw);
    }
    if (root->ne != NULL){
        freetree(root->ne);
    }
    if (root->sw != NULL){
        freetree(root->sw);
    }
    if (root->se != NULL){
        freetree(root->se);
    }

    free(root);

}

void mass_center(dummy_attribute* root, Node* level){
    if (level->quadrant_body_count == 1){
        root->mass += (*(level->bodies))->mass;
        root->position.x = (((*(level->bodies))->mass*(*(level->bodies))->position.x)+((root->mass-(*(level->bodies))->mass)*root->position.x))/root->mass;
        root->position.y = (((*(level->bodies))->mass*(*(level->bodies))->position.y)+((root->mass-(*(level->bodies))->mass)*root->position.y))/root->mass;

    }else{
        if(level->sw != NULL){
            mass_center(root,level->sw);
        }
        if(level->se != NULL){
            mass_center(root,level->se);
        }
        if(level->nw != NULL){
            mass_center(root,level->nw);
        }
        if(level->ne != NULL){
            mass_center(root,level->ne);
        }

    }


}

void create_bh(body_attribute* body, Node* root, dummy_attribute* bh_list,double theta_max){

    double dx =body->position.x-root->quadrant_center.x;
    double dy =body->position.y-root->quadrant_center.y;
    double distance = sqrt((dx*dx) + (dy*dy));
    double theta    = root->quadrant_size/distance;
    if(theta < theta_max || root->quadrant_body_count == 1){

        bh_list[bh_list[0].id].mass = root->quadrant_mass;
        bh_list[bh_list[0].id].position.x = root->quadrant_mass_center.x;
        bh_list[bh_list[0].id].position.y = root->quadrant_mass_center.y;
        bh_list[0].id++;

    }else{
        if (root->nw!=NULL){
            create_bh(body,root->nw,bh_list,theta_max);
        }
        if (root->ne!=NULL){
            create_bh(body,root->ne,bh_list,theta_max);
        }
        if (root->sw!=NULL){
            create_bh(body,root->sw,bh_list,theta_max);
        }
        if (root->se!=NULL){
            create_bh(body,root->se,bh_list,theta_max);
        }

    }
}///get this over with


void solver(body_attribute* body,dummy_attribute* cluster_list, double G, double eps,double dt){

    double fx = 0;
    double fy = 0;

    for (int i = 0; i < cluster_list[0].id; ++i) {

        double dx;
        double dy;
        double r;


        dx = cluster_list[i].position.x - body->position.x ;
        dy = cluster_list[i].position.y - body->position.y ;
        r = sqrt((dx*dx)+(dy*dy)) + eps;

        fx += G*(cluster_list[i].mass)*dx/(r*r*r);
        fy += G*(cluster_list[i].mass)*dy/(r*r*r);


    }

    body->velocity.x += dt*fx;
    body->velocity.y += dt*fy;

}


int main(int argc, char** argv) {

    ///arg check
    if(argc != 7){
        printf("wrong arguments\n");
        exit(1);
    }

    ///demangling
    int N            = atoi(argv[1]);
    char* filename   =     (argv[2]);
    int step_num     = atoi(argv[3]);
    double dt        = atof(argv[4]);
    double theta_max = atof(argv[5]);
    graphics     = atoi(argv[6]);


    ///filling in the bodies to the struct form
    FILE *fp = fopen(filename,"r");
    if (fp == NULL){
        printf("error in file opening\n");
        exit(1);
    }

    ///all bodies
    body_attribute* body_list[N];


    ///filling body list in
    int counter = 1;
    for(int i = 0; i<N;i++){
        body_list[i] = (body_attribute*) malloc(sizeof(body_attribute));
        for (int j = 0; j < 6; ++j) {
            double d = 0.0;
            fread(&d,1,sizeof(double),fp);
            if (counter == 1){
                (*body_list[i]).position.x = d;
            }else if(counter == 2){
                (*body_list[i]).position.y = d;
            }else if(counter == 3){
                (*body_list[i]).mass  = d;
            }else if(counter == 4){
                (*body_list[i]).velocity.x = d;
            }else if(counter == 5){
                (*body_list[i]).velocity.y = d;
            }else if(counter == 6){
                (*body_list[i]).brightness = d;
            }
            counter++;
        }
        counter = 1;
    }



    ///intialize graphics
    if (graphics){
        InitializeGraphics(argv[0],1100,1100);
        SetCAxes(0,1);
        ClearScreen();
        Refresh();
        ClearScreen();
        Refresh();
    }

    ///calculation of force
    double G  = 100.0/N;
    double eps = 0.001;






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int iteration = 0;
    while (iteration < step_num) {

        ///create the tree root
        Node* domain = NULL;
        domain = (Node*)calloc(1,sizeof(Node));
        domain->quadrant_center.x =  0.5;
        domain->quadrant_center.y =  0.5;
        domain->quadrant_size     =    1;
        domain->quadrant_body_count =  0;
        domain->quadrant_mass = 0;
        domain->quadrant_mass_center.x = 0;
        domain->quadrant_mass_center.y = 0;



        for (int i = 0; i < N; ++i) {
            ///create tree
            create_tree(&body_list[i],domain);
        }


        for (int ii = 0; ii < N; ++ii) {

            if(depth) DrawRectangle(domain->quadrant_center.x - (domain->quadrant_size)/2,
                                    domain->quadrant_center.y - (domain->quadrant_size)/2, W, H, domain->quadrant_size,
                                    domain->quadrant_size,0, yellowColor);
            if (graphics) {
                for (int i = 0; i < N; ++i) {
                    DrawCircle(body_list[i]->position.x, body_list[i]->position.y, W, H, 0.002, 0, redColor);
                }
            }

            ///initializing bh list
            dummy_attribute* bh_list = NULL;
            bh_list=(dummy_attribute*)calloc(N, sizeof(dummy_attribute));
            bh_list[0].id=0;
            create_bh(body_list[ii],domain, bh_list,theta_max);

            if (graphics){
                Refresh();
                usleep(6000);
                ClearScreen();
            }
            ///solver
            solver(body_list[ii],bh_list,G,eps,dt);

            free(bh_list);
        }
        ///calculate new locations body_list
        for (int i = 0; i < N; ++i) {

            (*body_list[i]).position.x = (*body_list[i]).position.x + (dt*(*body_list[i]).velocity.x);
            (*body_list[i]).position.y = (*body_list[i]).position.y + (dt*(*body_list[i]).velocity.y);

        }
        ///here
        ++iteration;
        freetree(domain);
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///write to file
    FILE* output;
    output = fopen("output.gal", "w+");
    for (int i = 0; i < N; ++i) {
        fwrite(&body_list[i]->position.x, 1, sizeof(double),output);
        fwrite(&body_list[i]->position.y, 1, sizeof(double),output);
        fwrite(&body_list[i]->mass,       1, sizeof(double),output);
        fwrite(&body_list[i]->velocity.x, 1, sizeof(double),output);
        fwrite(&body_list[i]->velocity.y, 1, sizeof(double),output);
        fwrite(&body_list[i]->brightness, 1, sizeof(double),output);


    }
    fclose(output);//as george sais, be a nice citizen of OS

    for (int i = 0; i < N; ++i) {
        free(body_list[i]);
    }

    return 0;
}///main