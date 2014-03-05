using System;
using System.Collections.Generic;
using System.Text;

namespace GeneticsLab
{
    class PairWiseAlign
    {
        
        /// <summary>
        /// Align only 5000 characters in each sequence.
        /// </summary>
        private int MaxCharactersToAlign = 5000; 

        /// <summary>
        /// this is the function you implement.
        /// </summary>
        /// <param name="sequenceA">the first sequence</param>
        /// <param name="sequenceB">the second sequence, may have length not equal to the length of the first seq.</param>
        /// <param name="resultTableSoFar">the table of alignment results that has been generated so far using pair-wise alignment</param>
        /// <param name="rowInTable">this particular alignment problem will occupy a cell in this row the result table.</param>
        /// <param name="columnInTable">this particular alignment will occupy a cell in this column of the result table.</param>
        /// <returns>the alignment score for sequenceA and sequenceB.  The calling function places the result in entry rowInTable,columnInTable
        /// of the ResultTable</returns>
        public int Align(GeneSequence sequenceA, GeneSequence sequenceB, ResultTable resultTableSoFar, int rowInTable, int columnInTable)
        {
            // Only fill in above the diagonal
            if (rowInTable <= columnInTable) return 0;

            string a = sequenceA.Sequence;
            string b = sequenceB.Sequence;

            int m = Math.Min(a.Length, 5000);
            int n = Math.Min(b.Length, 5000);
            int[,] E = new int[5001,5001];

            for (int i = 0; i <= m; i++)
            {
                E[i, 0] = 5 * i;
            }

            for (int j = 1; j <= n; j++)
            {
                E[0, j] = 5 * j;
            }

            for (int i = 1; i <= m; i++)
            {
                for (int j = 1; j <= n; j++)
                {
                    var indels = Math.Min(E[i - 1, j] + 5, E[i, j - 1] + 5);
                    var diff = (a[i - 1] == b[j - 1]) ? -3 : 1;
                    E[i, j] = Math.Min(indels, E[i - 1, j - 1] + diff);
                }
            }

            return E[m, n];
        }
    }
}
