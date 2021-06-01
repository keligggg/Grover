//Programmation de l'algorithme de recherche de Grover sur un problème de factorisation
//Science du numérique ISEP S2 05/2021
//Antoine Jeanneney, Adrien Sallé, Kelig Lefeuvre
//
//dotnet run --number _
//dotnet run --no-build --number _
//remplacer _ par les nombres

namespace GroversAlgorithm {
    //import des bibliothèques
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Preparation;

    @EntryPoint()
    operation factorisationDeGrover(number : Int) : Int {
        //Définition des paramètres
        //Definition de l'oracle
        let nbSolutions =4;
        let Oracle = chercheDiviseur(number, _, _);
        let phaseOracle = rebondDePhase(Oracle, _);
        //nombre de bit du nombre a factoriser
        let size = BitSizeI(number);
        //nombre d'itération optimal pour rapprocher la probabilité des 100% de bonne réponse
        let nbIterations = Round(PI() / 4.0 * Sqrt(IntAsDouble(size) / IntAsDouble(nbSolutions)));

        //initialisation du registre des Qubits
        use registre = Qubit[size] {
            iterationsGrover(registre, phaseOracle, nbIterations);
            let res = MultiM(registre);
            let answer = ResultArrayAsInt(res);
            Message($"The number {answer} is a factor of {number}.");
            //return pour l'analyse en python
            return answer;
        }

    }

    
    /// # Summary
    /// Exécute un nombre de fois donné l'algorithme de Grover
    /// # Input
    /// ## registre - Tableau de Qubit initialisés à l'état 0 
    /// 
    /// ## phaseOracle - Oracle de phase pour la tâche de Grover
    /// 
    /// ## iterations - Nombre d'itération correspondant à pi/4sqrt(nbElements/nbSolutions)
    /// 
    operation iterationsGrover(registre : Qubit[], phaseOracle : ((Qubit[]) => Unit is Adj), iterations : Int) : Unit {
        //Initialisation du registre
        ApplyToEach(H, registre);
        //Loop sur le nombre d'iterations
        for i in 1 .. iterations {
            //Application de l'oracle de phase
            phaseOracle(registre);
            //Application de l'opérateur de diffussion de Grover
            operateurDiffusion(registre);
        }
    }

    /// # Summary
    /// Opérateur de diffusion de Grover, permet d'inverser les états --> effet de miroir sur les amplitudes des Qubits
    /// # Input
    /// ## inputQubits - Tableau de Qubit
    /// 
    operation operateurDiffusion(inputQubits : Qubit[]) : Unit {
        within {
            ApplyToEachA(H, inputQubits);
            ApplyToEachA(X, inputQubits);
        } apply {
            Controlled Z(Most(inputQubits), Tail(inputQubits));
        }
    }


    //***Implémentation de l'oracle ***//

    /// # Summary
    /// Cherche un diviseur du nombre
    /// # Input
    /// ## dividende - Nombre de d'éléments
    /// 
    /// ## divisorregistre - Tableau de Qubit
    /// 
    /// ## target - Qubit cible qui donne le résultat, inversé si le nombre est un diviseur
    /// 
    operation chercheDiviseur (dividende : Int, divisorregistre : Qubit[], target : Qubit) : Unit is Adj+Ctl {
        //taille du dividende
        let size = BitSizeI(dividende);

        //Allocation de deux Qubits pour le dividende et le résultat
        use (dividendeQubits, resultQubits) = (Qubit[size], Qubit[size]);

        //Utilisation du Little Endian pour utiliser DivideI
        let reste = LittleEndian(dividendeQubits);
        let diviseur = LittleEndian(divisorregistre);
        let result = LittleEndian(resultQubits);

        within{
            //Encodage du dividende
            ApplyXorInPlace(dividende, reste);
            //Opération de division 
            DivideI(reste, diviseur, result);
            //Retourner tous les Qubits du reste
            ApplyToEachA(X, reste!);
        }
        apply{
            //Appliquez un Non controllé sur le reste et retourner l'état du Qubit cible si le reste est 0
            Controlled X(reste!, target);
        }
    }

    /// # Summary
    /// Algorithme du rebond de phase - Transforme un oracle de marquage en oracle de phase
    /// # Input
    /// ## Oracle - Oracle de phase
    /// 
    /// ## registre - Tableau de Qubit
    /// 
    operation rebondDePhase(Oracle : (Qubit[], Qubit) => Unit is Adj, registre : Qubit[]) : Unit is Adj {
        use target = Qubit();
        within {
            X(target);
            H(target);
        } apply {
            Oracle(registre, target);
        }
    }
}