$(document).ready(function() {
var tableName;

$('#update-button').on('click', function() {
    $('.date-div').show();
    tableName = $('input[name="dbOption"]:checked').val(); // db 선택
    var selectedFields = [];

    // 선택된 체크박스 필드 수집
    $('.data-checkbox:checked').each(function() {
        var $row = $(this).closest("tr"); // 현재 체크박스가 속한 행을 가져옴
        const data = $row.find("td[data-id]").data("id"); // 해당 행에서 'data-id' 값을 가져옴
        selectedFields.push(data); // 가져온 'data-id' 값을 필드로 사용
    });

    if (tableName) {
        $.getJSON('/get_columns/' + tableName + '/', function(data) {
            var $table = $('#data-table');
            $table.find('thead').empty(); // 테이블 헤더 초기화
            $table.find('tbody').empty(); // 테이블 바디 초기화

            // 필터링된 필드를 테이블에 추가
            var filteredFields = data.fields.filter(function(field) {
                return field.help_text.includes('check');
            });

            // 테이블 헤더 추가
            var $headerRow = $('<tr></tr>');
            $headerRow.append($('<th></th>').text('Chem_id'));
            $.each(filteredFields, function(index, field) {
                $headerRow.append($('<th></th>').text(field.name));
            });
            $table.find('thead').append($headerRow);

            // 사용자 입력을 위한 빈 행 추가
            $.each(selectedFields, function(index, fieldValue) {
                var $inputRow = $('<tr></tr>');
                $inputRow.append($('<td></td>').text(fieldValue));
                $.each(filteredFields, function(i, field) {
                    $inputRow.append($('<td></td>').append($('<input>').attr({
                        type: 'text',
                        name: field.name
                    })));
                });
                $table.find('tbody').append($inputRow);


            var defaultFields = data.fields.filter(function(field) {
                return field.help_text === 'default';
            });

            // 기존의 Comment div가 있으면 삭제
            $('.comment-div').remove();

            if (defaultFields.length > 0) {
                // 'default' help_text 가진 필드마다 Comment div 요소 추가

                $.each(defaultFields, function(index, field) {
                    var $commentDiv = $('<div></div>', {
                        class: 'row mb-3  comment-div'
                    }).append(
                        $('<label></label>', {
                            for: field.name,
                            class: 'col-sm-2 col-form-label',
                            text: field.name
                        })
                    ).append(
                        $('<div></div>', {
                            class: 'col-sm-10'
                        }).append(
                            $('<textarea></textarea>', {
                                class: 'form-control',
                                id:  field.name,
                                name: field.name,

                            })
                        )
                    );
                    // 추가한 div를 테이블 밑에 삽입
                    $table.before($commentDiv);
                     });
                } else {
                    console.log("No default field found.");
                }
            });
        }).fail(function(jqXHR, textStatus, errorThrown) {
            console.error('Error fetching columns:', textStatus, errorThrown);
        });
    } else {
        alert('Please select a table.');
        }
    });

$('#save-button').on('click', function() {
    var tableData = new FormData();
//  var tableData = [];
    var date = $('#date').val();
    var cell = $('#cell').val();
    var comment = $('#comment').val();
    var user = $('#id_userid').val();

    if( !date ){
    	alert("빈 칸을 입력해주세요.");
        return;
    }

   $('#data-table tbody tr').each(function(index, element) {  // index와 element를 추가
        var rowData = {};
        $(element).find('input').each(function() {
            var fieldName = $(this).attr('name');
            var fieldValue = $(this).val();
            rowData[fieldName] = fieldValue;
        });

        // 'selectedFields'에서 가져온 값들을 'rowData'에 추가
        var dataId = $(this).find('td').first().text(); // 첫 번째 셀의 텍스트를 'data-id'로 사용
        if (dataId) {
            rowData['chemical'] = dataId; // 'data_id' 추가
        }

        rowData['db_name'] = tableName; // 라디오박스 값을 추가
        rowData['date'] = date;
        rowData['cell'] = cell;
        rowData['comment'] = comment;
        rowData['user'] = user;

        tableData.append(`row${index}`, JSON.stringify(rowData));
        // index를 활용하여 key 지정
        //tableData.push(rowData);
    });
    // 이미지 파일 추가
    var images = $('#image')[0].files;
    $.each(images, function(index, file) {
        tableData.append('images', file);
    });

    $.ajax({
        type: 'POST',
        url: '/save_table_data/',
        data: tableData,
        processData: false,
        contentType: false,
        success: function(response) {
            alert('입력되었습니다.');
            location.reload();
        },
        error: function(error) {
            console.error('Error saving data:', error);
            alert('실패하였습니다.');
        }
    });
});
});
