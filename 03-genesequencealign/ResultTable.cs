using System;
using System.Collections.Generic;
using System.Text;
using System.Windows.Forms;

namespace GeneticsLab
{
    class ResultTable
    {
        DataGridView m_dataView;

        int m_numberOfSequences;

        double[,] m_results;

        public ResultTable(DataGridView dataView, int numberOfSequences)
        {
            m_dataView = dataView;
            m_numberOfSequences = numberOfSequences;

            m_results = new double[m_numberOfSequences, m_numberOfSequences];

            for (int i = 0; i < m_numberOfSequences; ++i)
            {
                DataGridViewTextBoxColumn column = new DataGridViewTextBoxColumn();
                column.MaxInputLength = 8;
                column.Width = 50;
                DataGridViewCell cell = new DataGridViewTextBoxCell();
                column.CellTemplate = cell;
                m_dataView.Columns.Add(column);
            }
            for (int j = 0; j < m_numberOfSequences; ++j)
            {
                DataGridViewRow row = new DataGridViewRow();
                row.Height = 20;
                for (int i = 0; i < m_numberOfSequences; ++i)
                {
                    DataGridViewCell cell = new DataGridViewTextBoxCell();
                    //cell.Value = ((double)i);
                    //cell.Value = "asdas";
                    row.Cells.Add(cell);
                }
                m_dataView.Rows.Add(row);
            }
        }

        public void SetCell(int x, int y, int value)
        {
            m_results[y, x] = value;
            m_dataView.Rows[y].Cells[x].Value = value;
            m_dataView.Update(); 
        }

        public void SetCell(int x, int y, string value)
        {
            m_results[y, x] = -1.0;
            m_dataView.Rows[y].Cells[x].Value = value;
            m_dataView.Update(); 
        }

        public double GetCell(int x, int y)
        {
            return (m_results[y, x]);
        }
    }
}
