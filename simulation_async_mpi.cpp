#include <cstdlib>
#include <random>
#include <iostream>
#include <fstream>
#include "contexte.hpp"
#include "individu.hpp"
#include "graphisme/src/SDL2/sdl2.hpp"
# include <chrono>
#include<mpi.h>
void màjStatistique( épidémie::Grille& grille, std::vector<épidémie::Individu> const& individus )
{
    for ( auto& statistique : grille.getStatistiques() )
    {
        statistique.nombre_contaminant_grippé_et_contaminé_par_agent = 0;
        statistique.nombre_contaminant_seulement_contaminé_par_agent = 0;
        statistique.nombre_contaminant_seulement_grippé              = 0;
    }
    auto [largeur,hauteur] = grille.dimension();
    auto& statistiques = grille.getStatistiques();
    for ( auto const& personne : individus )
    {
        auto pos = personne.position();

        std::size_t index = pos.x + pos.y * largeur;
        if (personne.aGrippeContagieuse() )
        {
            if (personne.aAgentPathogèneContagieux())
            {
                statistiques[index].nombre_contaminant_grippé_et_contaminé_par_agent += 1;
            }
            else 
            {
                statistiques[index].nombre_contaminant_seulement_grippé += 1;
            }
        }
        else
        {
            if (personne.aAgentPathogèneContagieux())
            {
                statistiques[index].nombre_contaminant_seulement_contaminé_par_agent += 1;
            }
        }
    }
}

void afficheSimulation(sdl2::window& écran, std::vector<épidémie::Grille::StatistiqueParCase> statt, std::size_t jour,int hauteur_grille,int largeur_grille)
//void afficheSimulation(sdl2::window& écran, épidémie::Grille const& grille, std::size_t jour)
{
    auto [largeur_écran,hauteur_écran] = écran.dimensions();
    
    //auto [largeur_grille,hauteur_grille] = grille.dimension();
    //auto const& statistiques = grille.getStatistiques();
    sdl2::font fonte_texte("./graphisme/src/data/Lato-Thin.ttf", 18);
    écran.cls({0x00,0x00,0x00});
    // Affichage de la grille :
    std::uint16_t stepX = largeur_écran/largeur_grille;
    unsigned short stepY = (hauteur_écran-50)/hauteur_grille;
    double factor = 255./15.;
    
    
    

    for ( unsigned short i = 0; i < largeur_grille; ++i )
    {
        for (unsigned short j = 0; j < hauteur_grille; ++j )
        {
            auto const& stat = statt[i+j*largeur_grille];
            int valueGrippe = stat.nombre_contaminant_grippé_et_contaminé_par_agent+stat.nombre_contaminant_seulement_grippé;
            int valueAgent  = stat.nombre_contaminant_grippé_et_contaminé_par_agent+stat.nombre_contaminant_seulement_contaminé_par_agent;
            std::uint16_t origx = i*stepX;
            std::uint16_t origy = j*stepY;
            std::uint8_t red = valueGrippe > 0 ? 127+std::uint8_t(std::min(128., 0.5*factor*valueGrippe)) : 0;
            std::uint8_t green = std::uint8_t(std::min(255., factor*valueAgent));
            std::uint8_t blue= std::uint8_t(std::min(255., factor*valueAgent ));
            écran << sdl2::rectangle({origx,origy}, {stepX,stepY}, {red, green,blue}, true);
        }
    }

    écran << sdl2::texte("Carte population grippée", fonte_texte, écran, {0xFF,0xFF,0xFF,0xFF}).at(largeur_écran/2, hauteur_écran-20);
    écran << sdl2::texte(std::string("Jour : ") + std::to_string(jour), fonte_texte, écran, {0xFF,0xFF,0xFF,0xFF}).at(0,hauteur_écran-20);
    écran << sdl2::flush;
    
  
}

void simulation(bool affiche,int rank,int nb,MPI_Comm newcomm)
{
	
MPI_Status status;

MPI_Request request;
    constexpr const unsigned int largeur_écran = 1280, hauteur_écran = 1024;
    sdl2::window écran("Simulation épidémie de grippe", {largeur_écran,hauteur_écran});



    if (rank != 0){
        sdl2::finalize();
    }

    unsigned int graine_aléatoire = 1;
    std::uniform_real_distribution<double> porteur_pathogène(0.,1.);


    épidémie::ContexteGlobal contexte;
    // contexte.déplacement_maximal = 1; <= Si on veut moins de brassage
    // contexte.taux_population = 400'000;
    //contexte.taux_population = 1'000;
    contexte.interactions.β = 60.;
    std::vector<épidémie::Individu> population;
    population.reserve(contexte.taux_population);
    épidémie::Grille grille{contexte.taux_population};

    auto [largeur_grille,hauteur_grille] = grille.dimension();
    // L'agent pathogène n'évolue pas et reste donc constant...
    épidémie::AgentPathogène agent(graine_aléatoire++);
    // Initialisation de la population initiale :
    for (std::size_t i = 0; i < contexte.taux_population; ++i )
    {
        std::default_random_engine motor(100*(i+1));
        population.emplace_back(graine_aléatoire++, contexte.espérance_de_vie, contexte.déplacement_maximal);
        population.back().setPosition(largeur_grille, hauteur_grille);
        if (porteur_pathogène(motor) < 0.2)
        {
            population.back().estContaminé(agent);   
        }
    }

    std::size_t jours_écoulés = 0;
    int         jour_apparition_grippe = 0;
    int         nombre_immunisés_grippe= (contexte.taux_population*23)/100;
    sdl2::event_queue queue;

    bool quitting = false;

    std::ofstream output("Courbe.dat");
    output << "# jours_écoulés \t nombreTotalContaminésGrippe \t nombreTotalContaminésAgentPathogène()" << std::endl;

    épidémie::Grippe grippe(0);


    std::cout << "Début boucle épidémie" << std::endl << std::flush;
    std::chrono::time_point < std::chrono::system_clock > start, end;


    
    
    while (!quitting)
    {
          
        if (rank==0){

        //start = std::chrono::system_clock::now();
      int flag=0;
        auto events = queue.pull_events();
        for ( const auto& e : events)
        {
            if (e->kind_of_event() == sdl2::event::quit)
                quitting = true;
        }

        //reçoi des infos
std::vector<épidémie::Grille::StatistiqueParCase> stat(largeur_grille*hauteur_grille);
//stat = grille.getStatistiques();


        
            MPI_Iprobe( 1, 0, MPI_COMM_WORLD, &flag, &status );
        if(flag){

std::cout<<"here"<<std::endl;
    
    MPI_Recv(&jours_écoulés , 1, MPI_INT,  1, 0, MPI_COMM_WORLD ,& status );
    
MPI_Recv(stat.data(), 3*largeur_grille*hauteur_grille, MPI_INT,  1, 0, MPI_COMM_WORLD ,&status );
        }
        //sinon ne5edh el 9aleb

if (affiche) 
       {
		   afficheSimulation(écran, stat, jours_écoulés,hauteur_grille,largeur_grille);
		}
        

//end = std::chrono::system_clock::now();//pour le calcul du temps d'affichage
//std::chrono::duration < double >elapsed_seconds = end - start;
//std::cout<<elapsed_seconds.count()<<std::endl;
        }


        


        else{
            //proc respo de la simutlation



        start = std::chrono::system_clock::now();

        if (jours_écoulés%365 == 0)// Si le premier Octobre (début de l'année pour l'épidémie ;-) )
        {
            grippe = épidémie::Grippe(jours_écoulés/365);
            jour_apparition_grippe = grippe.dateCalculImportationGrippe();
            grippe.calculNouveauTauxTransmission();
            // 23% des gens sont immunisés. On prend les 23% premiers
            for ( int ipersonne = 0; ipersonne < nombre_immunisés_grippe; ++ipersonne)
            {
                population[ipersonne].devientImmuniséGrippe();
            }
            for ( int ipersonne = nombre_immunisés_grippe; ipersonne < int(contexte.taux_population); ++ipersonne )
            {
                population[ipersonne].redevientSensibleGrippe();
            }
        }
        if (jours_écoulés%365 == std::size_t(jour_apparition_grippe))
        {
            for (int ipersonne = nombre_immunisés_grippe; ipersonne < nombre_immunisés_grippe + 25; ++ipersonne )
            {
                population[ipersonne].estContaminé(grippe);
            }
        }
        // Mise à jour des statistiques pour les cases de la grille :
        màjStatistique(grille, population);
        // On parcout la population pour voir qui est contaminé et qui ne l'est pas, d'abord pour la grippe puis pour l'agent pathogène
        std::size_t compteur_grippe = 0, compteur_agent = 0, mouru = 0;
//debut de la parallelisation des individus





        //for ( auto& personne : population )
int pas=population.size()/(nb-1);//on fait la simulation sut tous les procs sauf pour le proc 0

        for(int i=(rank-1)*pas;i<rank*pas;i++)//parallelisation de la boucle
        {

            if (population[i].testContaminationGrippe(grille, contexte.interactions, grippe, agent))
            {
                compteur_grippe ++;
                population[i].estContaminé(grippe);
            }
            if (population[i].testContaminationAgent(grille, agent))
            {
                compteur_agent ++;
                population[i].estContaminé(agent);
            }
            // On vérifie si il n'y a pas de personne qui veillissent de veillesse et on génère une nouvelle personne si c'est le cas.
            if (population[i].doitMourir())
            {
                mouru++;
                unsigned nouvelle_graine = jours_écoulés + population[i].position().x*population[i].position().y;
                population[i] = épidémie::Individu(nouvelle_graine, contexte.espérance_de_vie, contexte.déplacement_maximal);
                population[i].setPosition(largeur_grille, hauteur_grille);
            }
            population[i].veillirDUnJour();
            population[i].seDéplace(grille);
        }
//collecte des données


            MPI_Allgather(&population[(rank-1)*pas],pas*sizeof(épidémie::Individu),MPI_PACKED,&population[0],pas*sizeof(épidémie::Individu),MPI_PACKED,newcomm);
 int loccompteur_grippe = compteur_grippe;
            int loccompteur_agent = compteur_agent;
            int locmouru = mouru;
            MPI_Allreduce(&loccompteur_grippe, &compteur_grippe, 1, MPI_INT,MPI_SUM,newcomm);//j'ai utilise allreduce pas reduce car au cas on 2 procs en general çca marche pas
            MPI_Allreduce(&loccompteur_agent, &compteur_agent, 1, MPI_INT, MPI_SUM,newcomm);
            MPI_Allreduce(&locmouru, &mouru, 1, MPI_INT, MPI_SUM, newcomm);
        
if(rank==1) {
        //on va envoyer la grille et le jour a proc 0
        auto& stat = grille.getStatistiques();


        MPI_Isend(&jours_écoulés, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&request);
        MPI_Isend(stat.data(), 3*largeur_grille*hauteur_grille, MPI_INT, 0, 0, MPI_COMM_WORLD,&request);

}

        /*std::cout << jours_écoulés << "\t" << grille.nombreTotalContaminésGrippe() << "\t"
                  << grille.nombreTotalContaminésAgentPathogène() << std::endl;*/

        output << jours_écoulés << "\t" << grille.nombreTotalContaminésGrippe() << "\t"
               << grille.nombreTotalContaminésAgentPathogène() << std::endl;
        jours_écoulés += 1;

        end = std::chrono::system_clock::now();
std::chrono::duration < double >elapsed_seconds = end - start;
std::cout<<elapsed_seconds.count()<<std::endl;
        

        }


        
        
        
        
                
//std::chrono::duration < double >elapsed_seconds = end - start;
//std::cout<<elapsed_seconds.count()<<std::endl;








    }// Fin boucle temporelle
    output.close();
   
     
}



int main(int argc, char* argv[])
{


int nb , rank;
// Initialisation de MPI
MPI_Init ( &argc , &argv );
// Lit le nombre de tâches
MPI_Comm_size ( MPI_COMM_WORLD , & nb );
// Lit mon rang
MPI_Comm_rank ( MPI_COMM_WORLD , &rank);
// Ici on peut commencer à faire de la programmation
// parallèle par échange de messages .
// On termine l'environnement MPI


int color=rank>0;//pour isoler le proc dans un comm et les autres porc dans un autre comm

//Split the communicator based on the color and use the
// original rank for ordering
MPI_Comm newcomm;
MPI_Comm_split(MPI_COMM_WORLD, color, rank, &newcomm);




    // parse command-line
    bool affiche = true;
    for (int i=0; i<argc; i++) {
      std::cout << i << " " << argv[i] << "\n";
      if (std::string("-nw") == argv[i]) affiche = false;
    }
  
    //

    sdl2::init(argc, argv);
    {
		   
        simulation(affiche,rank,nb,newcomm);
        
    }
    
    sdl2::finalize();


MPI_Barrier(MPI_COMM_WORLD); 

MPI_Finalize ();

return MPI_SUCCESS ;
}
