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

        struct DpTableEntry
        {
            public int i;
            public int j;
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
            int[][] E = new int[MaxCharactersToAlign + 1][];
            E[0] = new int[MaxCharactersToAlign + 1];

            for (int j = 0; j <= n; j++)
            {
                E[0][j] = 5 * j;
            }

            for (int i = 1; i <= m; i++)
            {
                int active = i;
                int previous = i - 1;
                E[active] = new int[MaxCharactersToAlign + 1];
                for (int j = 0; j <= n; j++)
                {
                    var indels = E[previous][j] + 5;

                    if (j == 0)
                    {
                        E[active][j] = indels;
                    }
                    else
                    {
                        indels = Math.Min(E[active][j - 1] + 5, indels);
                        var diff = (a[i - 1] == b[j - 1]) ? -3 : 1;
                        E[active][j] = Math.Min(indels, E[previous][j - 1] + diff);
                    }


                }
            }

            return E[m][n];
        }

        public string[] Extract(GeneSequence sequenceA, GeneSequence sequenceB)
        {
            string a = sequenceA.Sequence;
            string b = sequenceB.Sequence;
            int m = Math.Min(a.Length, MaxCharactersToExtract);
            int n = Math.Min(b.Length, MaxCharactersToExtract);
            int[,] E = new int[MaxCharactersToExtract + 1, MaxCharactersToExtract + 1];
            Dictionary<DpTableEntry, Direction> prev = new Dictionary<DpTableEntry, Direction>();

            // (0, 0) has 0 edit distance and no "previous" pointer
            E[0, 0] = 0;
            prev.Add(new DpTableEntry { i = 0, j = 0 }, Direction.NONE);

            for (int i = 1; i <= m; i++)
            {
                E[i, 0] = 5 * i;
                prev.Add(new DpTableEntry { i = i, j = 0 }, Direction.VERTICAL);
            }

            for (int j = 1; j <= n; j++)
            {
                E[0, j] = 5 * j;
                prev.Add(new DpTableEntry { i = 0, j = j }, Direction.HORIZONTAL);
            }

            for (int i = 1; i <= m; i++)
            {
                for (int j = 1; j <= n; j++)
                {
                    var vertical = E[i - 1, j] + 5;
                    var horizontal = E[i, j - 1] + 5;
                    var diagonal = E[i - 1, j - 1] + ((a[i - 1] == b[j - 1]) ? -3 : 1);
                    if (diagonal <= vertical && diagonal <= horizontal)
                    {
                        E[i, j] = diagonal;
                        prev.Add(new DpTableEntry { i = i, j = j }, Direction.DIAGONAL);
                    }
                    else if (vertical <= horizontal)
                    {
                        E[i, j] = vertical;
                        prev.Add(new DpTableEntry { i = i, j = j }, Direction.VERTICAL);
                    }
                    else
                    {
                        E[i, j] = horizontal;
                        prev.Add(new DpTableEntry { i = i, j = j }, Direction.HORIZONTAL);
                    }
                }
            }

            StringBuilder alignmentA = new StringBuilder();
            StringBuilder alignmentB = new StringBuilder();
            DpTableEntry current = new DpTableEntry { i = m, j = n };
            Direction next;
            while ((next = prev[current]) != Direction.NONE)
            {
                if (next == Direction.DIAGONAL)
                {
                    current = new DpTableEntry { i = current.i - 1, j = current.j - 1 };
                    alignmentA.Append(a[current.i]);
                    alignmentB.Append(b[current.j]);
                }
                else if (next == Direction.VERTICAL)
                {
                    current = new DpTableEntry { i = current.i - 1, j = current.j };
                    alignmentA.Append(a[current.i]);
                    alignmentB.Append('-');
                }
                else
                {
                    current = new DpTableEntry { i = current.i, j = current.j - 1 };
                    alignmentA.Append('-');
                    alignmentB.Append(b[current.j]);
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
