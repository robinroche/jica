/*	
 * Copyright 2011, Robin Roche
 * This file is part of jica.

    jica is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    jica is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with jica.  If not, see <http://www.gnu.org/licenses/>.
*/



/**
 * A simple test class for jica
 * @author Robin Roche
 */
public class TestICA 
{

	/**
	 * @param args
	 */
	public static void main(String[] args) 
	{
		int problemDimension = 2; 
		double[] minBounds = {-10,-10};
		double[] maxBounds = {10,10};
		
		Object[] argsICA ={problemDimension,minBounds,maxBounds};
		ICAlgorithm ica = new ICAlgorithm(argsICA);

		ica.runICA();
		
	}

}
