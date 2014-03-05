using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Text;
using System.Windows.Forms;
using System.Diagnostics;


namespace GeneticsLab
{
    public partial class MainForm : Form
    {
        DatabaseController m_dbController;

        ResultTable m_resultTable;

        GeneSequence[] m_sequences;

        public MainForm()
        {
            InitializeComponent();

            m_dbController = new DatabaseController();
            m_dbController.EstablishConnection("../../db1.mdb");

            statusMessage.Text = "Loading Database...";

            // Set the number of Sequences to load below.
            m_sequences = m_dbController.ReadGeneSequences(10);

            m_resultTable = new ResultTable(this.dataGridViewResults, m_sequences.Length);

            statusMessage.Text = "Loaded Database.";

        }

        private void fillMatrix()
        {
            PairWiseAlign processor = new PairWiseAlign();
            for (int y = 0; y < m_sequences.Length; ++y)
            {
                for (int x = 0; x < m_sequences.Length; ++x)
                {
                    m_resultTable.SetCell(x, y, processor.Align(m_sequences[x], m_sequences[y],m_resultTable,x,y));
                    //m_resultTable.SetCell(x, y, ("(" + x + ", " + y + ")"));
                }
            }
        }

        private void processButton_Click(object sender, EventArgs e)
        {
            statusMessage.Text = "Processing...";
            Stopwatch timer = new Stopwatch();
            timer.Start();
                   fillMatrix();
            timer.Stop();
            statusMessage.Text = "Done.  Time taken: " + timer.Elapsed;

        }
    }
}