SELECT exists(
  SELECT *

  FROM `document_annotation`

  INNER JOIN `relation_relationannotation`
    ON `relation_relationannotation`.`id` = `document_annotation`.`object_id` AND `relation_relationannotation`.`relation_id` = {relation_id}

  INNER JOIN `document_view`
    ON `document_view`.`id` = `document_annotation`.`view_id` AND `document_view`.`user_id` = {user_id}

  INNER JOIN `document_section`
    ON `document_section`.`id` = `document_view`.`section_id` AND `document_section`.`document_id` = {document_id}

  WHERE `document_annotation`.`kind` = 'r'
    AND `document_annotation`.`content_type_id` = 56
) as `exists`
