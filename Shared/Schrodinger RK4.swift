//
//  Schrodinger RK4.swift
//  Test Plot Threaded
//
//  Created by Matthew Malaker on 2/25/22.
//

import Foundation
import SwiftUI
import CorePlot
import Accelerate


class Schrodinger: NSObject, ObservableObject {
 
    var potential = Potential()
//    var potentialV: (onedArray: [Double], xArray: [Double], yArray: [Double]) = ([],[],[])
    
    
    var hbar2over2m = 0.0
    
    //k and l are the same type of variable but are used in the separate RK4 solutions
    //we need 2 solutions because we have a second order equation and need to solve for psi' and psi, not just psi
    //We do this because we only have a relation for psi' from psi'', which itself is a variation of psi, but psi' is unknown
    //past the initial condition.
//    var k1: Double = 0.0
//    var k2: Double = 0.0
//    var k3: Double = 0.0
//    var k4: Double = 0.0
//    var l1: Double = 0.0
//    var l2: Double = 0.0
//    var l3: Double = 0.0
//    var l4: Double = 0.0
    
    
//      2        2
//- hbar  partial psi
//-------- ----------- + VPsi = EPsi
//   2m              2
//          partial x

    

    
    
//      2        2
//- hbar  partial psi
//-------- ----------- + (V-E)Psi = 0
//   2m              2
//          partial x

    
    
// We need to solve for E, Psi, and the functional
    
//
//Psi'' = (2m(V-E)/hbar^2)*Psi
//
//
//
//
    
    override init() {
        
        
        //Must call super init before initializing plot
        super.init()
        
        let hbar = 6.582119569e-16
        let hbar2 = pow((6.582119569e-16),2.0)
        
        let massE = 510998.946/(pow(299792458.0e10,2.0))
        
        hbar2over2m = hbar2/(2.0*massE)
//      print(hbar2over2m)
    }

    func setV(choice: String, xMin: Double, xMax: Double, stepSize: Double){
        potential.getPotential(potentialToGet: choice, xMin: xMin, xMax: xMax, stepSize: stepSize)
        
    }
    

    func matrixSolve(xMin: Double, xMax: Double, stepSize: Double, potentialStr: String, E: Double, N: Int) async->[(eValue: (real: Double, imag: Double), eVector:[(real: Double, imag: Double)])]{
        potential.getPotential(potentialToGet: potentialStr, xMin: xMin, xMax: xMax, stepSize: stepSize)
//        var potentialV = Potential()
//        potentialV.getPotential(potentialToGet: potentialStr, xMin: xMin, xMax: xMax, stepSize: stepSize)
        
        
        print("stepSize = \(stepSize)")
        
        
        var H: [[Double]] = [[Double]](repeating: [Double](repeating: 0.0, count: N), count: N)
        let maxIndex = Int((xMax-xMin)/stepSize)
        var firstTerm: Double = 0.0
        var secondTerm = 0.0
        var valueForSecondTermIntegration = 0.0
        print(N)
        print(H.count)
        print(H[0].count)
        let a = xMax-xMin
        var arrayOfEigens: [(eValue: (real: Double, imag: Double), eVector:[(real: Double, imag: Double)])] = []
//        print(potentialV.potential.yArray)
        
        //Constructing H for each Hij
        //Threading this makes it a lot faster, and each (i,j) can be threaded
        
//        let HResults = await withTaskGroup(of: (i: Int, j: Int, Hij: Double).self, returning: [(i: Int, j: Int, Hij: Double)].self, body: {taskGroup in
//
//            for i in stride(from: 0, to: N, by: 1){
//
//                for j in stride(from: 0, to: N, by: 1){
//
//                    taskGroup.addTask{
//                        var firstTerm: Double = 0.0
//                        var secondTerm = 0.0
//                        var valueForSecondTermIntegration = 0.0
//                    if(i==j){
//                        firstTerm = pow(Double(j+1)*Double.pi/(xMax-xMin),2)*self.hbar2over2m // jth energy if i==j; 0 else
//                    }
//
//
//                    for index in stride(from: 0, to: maxIndex, by: 1){
//                        let x = Double(index)*stepSize + xMin
//
//                        valueForSecondTermIntegration += (2/(xMax-xMin))*sin(Double(i+1)*Double.pi*x/a)*self.potential.potential.yArray[index]*sin(Double(j+1)*Double.pi*x/a)
//
//
//                    }
//                            secondTerm = valueForSecondTermIntegration/(Double(maxIndex)/stepSize)
//
//
//                        return (i: i, j: j, Hij: firstTerm+secondTerm)
//                    firstTerm = 0.0
//                    secondTerm = 0.0
//                    }
//                }
//            }
//
//            var interimResults = [(i: Int, j: Int, Hij: Double)]()
//            //reordering results as they come in is not necessary
//            for await result in taskGroup{
//                interimResults.append(result)
//            }
//
//            return interimResults
//        })
//
//        for element in HResults{
//            H[element.i][element.j] = element.Hij
//        }
        
        print("maxIndex = \(maxIndex)")
        for i in stride(from: 0, to: N, by: 1){

            for j in stride(from: 0, to: N, by: 1){

                if(i==j){
                    firstTerm = pow(Double(j+1)*Double.pi/(xMax-xMin),2)*hbar2over2m // jth energy if i==j; 0 else
                }
                
                
                for index in stride(from: 0, to: maxIndex, by: 1){
                    let x = Double(index)*stepSize + xMin

                    valueForSecondTermIntegration += (2/(xMax-xMin))*sin(Double(i+1)*Double.pi*x/a)*potential.potential.yArray[index]*sin(Double(j+1)*Double.pi*x/a)

                    
                }
                        secondTerm = valueForSecondTermIntegration/(Double(maxIndex)/stepSize)


                H[i][j] = (firstTerm+secondTerm)
                firstTerm = 0.0
                secondTerm = 0.0
                
            }
        }
//        print(H)
        
        
        
        //Solve EigenValue problem for Hij
        arrayOfEigens = calculateEigenvalues(arrayForDiagonalization: pack2dArray(arr: H, rows: N, cols: N))
        arrayOfEigens = arrayOfEigens.sorted(by: {$0.eValue.real < $1.eValue.real})
        return arrayOfEigens
        
    }


    

    func normalize(data: (E: Double, psi: [Double]), stepSize: Double) -> (E: Double, psi: [Double]){
        var psiarray = data.psi
        var normarray: [Double] = []
        var normed: (E: Double, psi: [Double]) = (E: data.E, [])
        var integral: Double = 0.0
        
        //We need to integrate the two wavefunctions
        
        for i in psiarray{
            integral += pow(i,2.0)*stepSize
        }
        for i in psiarray{
            normarray.append(i*pow(integral,-0.5))
        }
        
        normed = (data.E, normarray)
        return normed
    }
    
    
    ///Function Courtesy Dr. Jeffery Terry
    /// pack2DArray
    /// Converts a 2D array into a linear array in FORTRAN Column Major Format
    ///
    /// - Parameters:
    ///   - arr: 2D array
    ///   - rows: Number of Rows
    ///   - cols: Number of Columns
    /// - Returns: Column Major Linear Array
    func pack2dArray(arr: [[Double]], rows: Int, cols: Int) -> [Double] {
        var resultArray = Array(repeating: 0.0, count: rows*cols)
        for Iy in 0...cols-1 {
            for Ix in 0...rows-1 {
                let index = Iy * rows + Ix
                resultArray[index] = arr[Ix][Iy]
            }
        }
        return resultArray
    }
    
    /// unpack2DArray
    /// Converts a linear array in FORTRAN Column Major Format to a 2D array in Row Major Format
    ///
    /// - Parameters:
    ///   - arr: Column Major Linear Array
    ///   - rows: Number of Rows
    ///   - cols: Number of Columns
    /// - Returns: 2D array
    func unpack2dArray(arr: [Double], rows: Int, cols: Int) -> [[Double]] {
        var resultArray = [[Double]](repeating:[Double](repeating:0.0 ,count:rows), count:cols)
        for Iy in 0...cols-1 {
            for Ix in 0...rows-1 {
                let index = Iy * rows + Ix
                resultArray[Ix][Iy] = arr[index]
            }
        }
        return resultArray
    }
    
    ///Function Courtesy Dr. Jeffery Terry, modified for use here
    /// calculateEigenvalues
    ///
    /// - Parameter arrayForDiagonalization: linear Column Major FORTRAN Array for Diagonalization
    /// - Returns: String consisting of the Eigenvalues and Eigenvectors
    func calculateEigenvalues(arrayForDiagonalization: [Double]) -> [(eValue:(real:Double, imag: Double), eVector:[(real: Double, imag: Double)])] {
        /* Integers sent to the FORTRAN routines must be type Int32 instead of Int */
        //var N = Int32(sqrt(Double(startingArray.count)))
        
        var returnString = ""
        var eigenVectorAtIndex: [(real: Double, imag: Double)] = []
        var returnArray: [(eValue:(real:Double, imag: Double), eVector:[(real: Double, imag: Double)])] = []
        
        var N = Int32(sqrt(Double(arrayForDiagonalization.count)))
        var N2 = Int32(sqrt(Double(arrayForDiagonalization.count)))
        var N3 = Int32(sqrt(Double(arrayForDiagonalization.count)))
        var N4 = Int32(sqrt(Double(arrayForDiagonalization.count)))
        
        var flatArray = arrayForDiagonalization
        
        var error : Int32 = 0
        var lwork = Int32(-1)
        // Real parts of eigenvalues
        var wr = [Double](repeating: 0.0, count: Int(N))
        // Imaginary parts of eigenvalues
        var wi = [Double](repeating: 0.0, count: Int(N))
        // Left eigenvectors
        var vl = [Double](repeating: 0.0, count: Int(N*N))
        // Right eigenvectors
        var vr = [Double](repeating: 0.0, count: Int(N*N))
        
        
        /* Eigenvalue Calculation Uses dgeev */
        /*   int dgeev_(char *jobvl, char *jobvr, Int32 *n, Double * a, Int32 *lda, Double *wr, Double *wi, Double *vl,
         Int32 *ldvl, Double *vr, Int32 *ldvr, Double *work, Int32 *lwork, Int32 *info);*/
        
        /* dgeev_(&calculateLeftEigenvectors, &calculateRightEigenvectors, &c1, AT, &c1, WR, WI, VL, &dummySize, VR, &c2, LWork, &lworkSize, &ok)    */
        /* parameters in the order as they appear in the function call: */
        /* order of matrix A, number of right hand sides (b), matrix A, */
        /* leading dimension of A, array records pivoting, */
        /* result vector b on entry, x on exit, leading dimension of b */
        /* return value =0 for success*/
        
        
        
        /* Calculate size of workspace needed for the calculation */
        
        var workspaceQuery: Double = 0.0
        dgeev_(UnsafeMutablePointer(mutating: ("N" as NSString).utf8String), UnsafeMutablePointer(mutating: ("V" as NSString).utf8String), &N, &flatArray, &N2, &wr, &wi, &vl, &N3, &vr, &N4, &workspaceQuery, &lwork, &error)
        
        print("Workspace Query \(workspaceQuery)")
        
        /* size workspace per the results of the query */
        
        var workspace = [Double](repeating: 0.0, count: Int(workspaceQuery))
        lwork = Int32(workspaceQuery)
        
        /* Calculate the size of the workspace */
        
        dgeev_(UnsafeMutablePointer(mutating: ("N" as NSString).utf8String), UnsafeMutablePointer(mutating: ("V" as NSString).utf8String), &N, &flatArray, &N2, &wr, &wi, &vl, &N3, &vr, &N4, &workspace, &lwork, &error)
        
        //The original function returns a string, which is not usable in this code. We need to have a structure containing the eigenvectors and eigenvalues for the potential. We will probably have an array of tuples
        if (error == 0)
        {
            for index in 0..<wi.count      /* transform the returned matrices to eigenvalues and eigenvectors */
            {
                
                for j in 0..<N
                {
                    
                    
                    if(wi[index]==0)
                    {
                        
                        //returnString += "\(vr[Int(index)*(Int(N))+Int(j)]) + 0.0i, \n" /* print x */
                        eigenVectorAtIndex.append((real: Double(vr[Int(index)*(Int(N))+Int(j)]), imag: Double(0.0)))
                        
                    }
                    else if(wi[index]>0)
                    {
                        if(vr[Int(index)*(Int(N))+Int(j)+Int(N)]>=0)
                        {
                            //returnString += "\(vr[Int(index)*(Int(N))+Int(j)]) + \(vr[Int(index)*(Int(N))+Int(j)+Int(N)])i, \n"
                            eigenVectorAtIndex.append((real: Double(vr[Int(index)*(Int(N))+Int(j)]), imag: Double(vr[Int(index)*(Int(N))+Int(j)+Int(N)])))
                        }
                        else
                        {
                            //returnString += "\(vr[Int(index)*(Int(N))+Int(j)]) - \(fabs(vr[Int(index)*(Int(N))+Int(j)+Int(N)]))i, \n"
                            eigenVectorAtIndex.append((real: Double(vr[Int(index)*(Int(N))+Int(j)]), imag: Double(vr[Int(index)*(Int(N))+Int(j)+Int(N)])))
                            
                        }
                    }
                    else
                    {
                        if(vr[Int(index)*(Int(N))+Int(j)]>0)
                        {
                            //returnString += "\(vr[Int(index)*(Int(N))+Int(j)-Int(N)]) - \(vr[Int(index)*(Int(N))+Int(j)])i, \n"
                            eigenVectorAtIndex.append((real: Double(vr[Int(index)*(Int(N))+Int(j)-Int(N)]), imag: Double(vr[Int(index)*(Int(N))+Int(j)])))
                        }
                        else
                        {
                            //returnString += "\(vr[Int(index)*(Int(N))+Int(j)-Int(N)]) + \(fabs(vr[Int(index)*(Int(N))+Int(j)]))i, \n"
                            eigenVectorAtIndex.append((real: Double(vr[Int(index)*(Int(N))+Int(j)-Int(N)]), imag: Double(fabs(vr[Int(index)*(Int(N))+Int(j)]))))
                            
                        }
                    }
                }
                
                

                    //returnString += "Eigenvalue\n\(wr[index]) + \(wi[index])i\n\n"
                    returnArray.append((eValue: (real:Double(wr[index]), imag:Double(wi[index])), eVector: eigenVectorAtIndex))
                    eigenVectorAtIndex.removeAll()
                
                //returnString += "Eigenvector\n"
                //returnString += "["
                
                
                /* To Save Memory dgeev returns a packed array if complex */
                /* Must Unpack Properly to Get Correct Result
                 
                 VR is DOUBLE PRECISION array, dimension (LDVR,N)
                 If JOBVR = 'V', the right eigenvectors v(j) are stored one
                 after another in the columns of VR, in the same order
                 as their eigenvalues.
                 If JOBVR = 'N', VR is not referenced.
                 If the j-th eigenvalue is real, then v(j) = VR(:,j),
                 the j-th column of VR.
                 If the j-th and (j+1)-st eigenvalues form a complex
                 conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
                 v(j+1) = VR(:,j) - i*VR(:,j+1). */
                
                
                
                /* Remove the last , in the returned Eigenvector */
                //returnString.remove(at: //returnString.index(before: //returnString.endIndex))
                //returnString.remove(at: //returnString.index(before: //returnString.endIndex))
                //returnString.remove(at: //returnString.index(before: //returnString.endIndex))
                //returnString += "]\n\n"
            }
        }
        else {print("An error occurred\n")}
        
        return (returnArray)
    }
    
    
    
    
    
}
