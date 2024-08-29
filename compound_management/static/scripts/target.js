$(document).ready(function() {
    $('#select_all').click(function() {
        $('input[name="selected_chems"]').prop('checked', this.checked);
    });

    const isChecked = localStorage.getItem('checkboxStatus') === 'true';
    $('#myCheckbox').prop('checked', isChecked);

    $('#delete_selected').click(function() {
        var selected = [];
        $('input[name="selected_chems"]:checked').each(function() {
            selected.push($(this).val());
        });

        if (selected.length > 0) {
            if (confirm('Are you sure you want to delete the selected chemicals?')) {
                $.ajax({
                    url: "{% url 'delete_selected_chems' target=target %}",
                    method: "POST",
                    data: {
                        'selected_chems': selected,
                        'csrfmiddlewaretoken': '{{ csrf_token }}'
                    },
                    success: function(response) {
                        location.reload();
                    },
                    error: function(xhr, status, error) {
                        alert('An error occurred while deleting the selected chemicals: ' + error);
                    }
                });
            }
        } else {
            alert('Please select at least one chemical to delete.');
        }
    });
    $('#view_selected').click(function() {
        var selected = [];
        $('input[name="selected_chems"]:checked').each(function() {
            selected.push($(this).val());
        });
         console.log(selected);
    });

    $('#select_check').on('change', function() {
        localStorage.setItem('checkboxStatus', $(this).is(':checked'));
        });

    // 이미지 클릭 이벤트 처리
    $(document).on('click', '.image-popup', function(event) {
        var row = $(this).closest('tr');
        var imgWindow = window.open("", "Image", "width=800,height=700, scrollbars=yes");

        imgWindow.document.write(`
            <!DOCTYPE html>
            <html lang="en">
            <head>
                <meta charset="UTF-8">
                <meta name="viewport" content="width=device-width, initial-scale=1.0">
                <title>Image</title>
                <style>
                     body {
                        margin: 0;
                        background-color: #f5f5f5;
                        display: flex;
                        align-items: center;
                        padding: 20px;
                    }
                    .container {
                        max-width: 800px;
                        width: 100%;
                        background: #fff;
                        box-shadow: 0 4px 8px rgba(0,0,0,0.1);
                        border-radius: 8px;
                        overflow: hidden;
                        padding: 20px;
                        box-sizing: border-box;
                        overflow-x: auto;
                        min-width: 30%; /* 콘텐츠가 최소한 화면 너비를 차지하도록 설정 */
                    }
                    img {
                        max-width: 100%;
                        height: auto;
                        display: block;
                        margin: 0 auto;
                        border-radius: 8px;
                    }
                    table {
                        width: 100%;
                        border-collapse: collapse;
                        margin-top: 20px;
                    }
                    th, td {
                        border: 1px solid #ddd;
                        padding: 8px;
                        text-align: left;
                    }
                    th {
                        background-color: #f4f4f4;
                    }
                    h1 {
                        margin-top: 0;
                    }

                </style>
                </head>
                <div class="container">
                    <img src="${this.src}" alt="${this.alt}">
                    <h1>{{ target }} / Chem_id : ${row.find('td').eq(2).text()}</h1>
                    <table>
                        <tr>
                            <td>MW</td>
                            <td>${row.find('td').eq(4).text()}</td>
                        </tr>
                        <tr>
                            <td>cLogP</td>
                            <td>${row.find('td').eq(5).text()}</td>
                        </tr>
                        <tr>
                            <td>TPSA</td>
                            <td>${row.find('td').eq(6).text()}</td>
                        </tr>
                        <tr>
                            <td>H_donors</td>
                            <td>${row.find('td').eq(7).text()}</td>
                        </tr>
                         <tr>
                            <td>H_acceptors</td>
                            <td>${row.find('td').eq(8).text()}</td>
                        </tr>
                         <tr>
                            <td>lipinski</td>
                            <td>${row.find('td').eq(9).text()}</td>
                        </tr>
                    </table>
                </div>
            </html>
        `);
    });
});

$(".favorite-btn").on("click", function() {
    var $button = $(this);
    var chemicalId = $button.data("id");
    var csrfToken = $('input[name=csrfmiddlewaretoken]').val();
    var isFavorite = $button.text().trim() === '❤️';

    $.ajax({
        url: `/item/${chem_id}/toggle_favorite/`,
        type: 'POST',
        contentType: 'application/json',
        headers: {
            'X-CSRFToken': csrfToken
        },
        data: JSON.stringify({ is_favorite: isFavorite }),
        success: function(data) {
            if (data.success) {
                $button.text(data.is_favorite ? '❤️' : '♡');
            } else {
                alert("An error occurred. Please try again.");
            }
        },
        error: function() {
            alert("An error occurred. Please try again.");
        }
    });
});


function getSortOrder(currentOrder, column) {
    if (currentOrder === 'asc' || currentOrder === '') {
        return 'desc';
    } else {
        return 'asc';
    }
}