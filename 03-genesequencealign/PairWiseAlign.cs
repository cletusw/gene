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
        private const int MaxCharactersToAlign = 5000;

        private const int MaxCharactersToExtract = 100;

        enum Direction
        {
            VERTICAL,
            HORIZONTAL,
            DIAGONAL,
            NONE
        }

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
            int m = Math.Min(a.Length, MaxCharactersToAlign);
            int n = Math.Min(b.Length, MaxCharactersToAlign);
            int[][] E = new int[2][];
            E[0] = new int[MaxCharactersToAlign + 1];
            E[1] = new int[MaxCharactersToAlign + 1];

            // Initialize first row with cost of indels
            for (int j = 0; j <= n; j++)
            {
                E[0][j] = 5 * j;
            }

            int previous = 0;
            int active = 1;

            for (int i = 1; i <= m; i++)
            {
                for (int j = 0; j <= n; j++)
                {
                    var indels = E[previous][j] + 5;

                    if (j == 0)
                    {
                        // If first element, only have one option
                        E[active][j] = indels;
                    }
                    else
                    {
                        // Get the minimum cost from the two available indels and the diagonal match/sub
                        indels = Math.Min(E[active][j - 1] + 5, indels);
                        var diff = (a[i - 1] == b[j - 1]) ? -3 : 1;
                        E[active][j] = Math.Min(indels, E[previous][j - 1] + diff);
                    }
                }

                // Swap active and previous (new active will get overwritten)
                active = (active == 0) ? 1 : 0;
                previous = (previous == 0) ? 1 : 0;
            }

            // Return last element in last filled row
            return E[previous][n];
        }

        public string[] Extract(GeneSequence sequenceA, GeneSequence sequenceB)
        {
            string a = sequenceA.Sequence;
            string b = sequenceB.Sequence;
            int m = Math.Min(a.Length, MaxCharactersToExtract);
            int n = Math.Min(b.Length, MaxCharactersToExtract);
            int[,] E = new int[MaxCharactersToExtract + 1, MaxCharactersToExtract + 1];
            Dictionary<Tuple<int, int>, Direction> prev = new Dictionary<Tuple<int, int>, Direction>();

            // (0, 0) has 0 edit distance and no "previous" pointer
            E[0, 0] = 0;
            prev.Add(Tuple.Create(0, 0), Direction.NONE);

            // Initialize first column with cost of indels
            for (int i = 1; i <= m; i++)
            {
                E[i, 0] = 5 * i;
                prev.Add(Tuple.Create(i, 0), Direction.VERTICAL);
            }

            // Initialize first row with cost of indels
            for (int j = 1; j <= n; j++)
            {
                E[0, j] = 5 * j;
                prev.Add(Tuple.Create(0, j), Direction.HORIZONTAL);
            }

            for (int i = 1; i <= m; i++)
            {
                for (int j = 1; j <= n; j++)
                {
                    var vertical = E[i - 1, j] + 5;
                    var horizontal = E[i, j - 1] + 5;
                    var diagonal = E[i - 1, j - 1] + ((a[i - 1] == b[j - 1]) ? -3 : 1);

                    // Get the minimum cost from the two available indels and the diagonal match/sub
                    if (diagonal <= vertical && diagonal <= horizontal)
                    {
                        E[i, j] = diagonal;
                        prev.Add(Tuple.Create(i, j), Direction.DIAGONAL);
                    }
                    else if (vertical <= horizontal)
                    {
                        E[i, j] = vertical;
                        prev.Add(Tuple.Create(i, j), Direction.VERTICAL);
                    }
                    else
                    {
                        E[i, j] = horizontal;
                        prev.Add(Tuple.Create(i, j), Direction.HORIZONTAL);
                    }
                }
            }

            // Store the alignment by walking backward from the final entry using prev pointers
            StringBuilder alignmentA = new StringBuilder();
            StringBuilder alignmentB = new StringBuilder();
            var current = Tuple.Create(m, n);
            Direction next;
            while ((next = prev[current]) != Direction.NONE)
            {
                if (next == Direction.DIAGONAL)
                {
                    current = Tuple.Create(current.Item1 - 1, current.Item2 - 1);

                    // Match or Substitution
                    alignmentA.Append(a[current.Item1]);
                    alignmentB.Append(b[current.Item2]);
                }
                else if (next == Direction.VERTICAL)
                {
                    current = Tuple.Create(current.Item1 - 1, current.Item2);

                    // Delete from a
                    alignmentA.Append(a[current.Item1]);
                    alignmentB.Append('-');
                }
                else
                {
                    current = Tuple.Create(current.Item1, current.Item2 - 1);

                    // Insert into b
                    alignmentA.Append('-');
                    alignmentB.Append(b[current.Item2]);
                }
            }

            // Reverse the alignment since we built it backward
            char[] alignmentAChar = alignmentA.ToString().ToCharArray();
            char[] alignmentBChar = alignmentB.ToString().ToCharArray();
            Array.Reverse(alignmentAChar);
            Array.Reverse(alignmentBChar);
            string[] alignment = new string[2];
            alignment[0] = new string(alignmentAChar);
            alignment[1] = new string(alignmentBChar);

            return alignment;
        }
    }
}
